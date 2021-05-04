###############################################
# Projections using transmission model and geostatistical map of Trachoma
#
#  Using algorithm written by Retkute et al. 2020
#  "Integrating geostatistical maps and transmission models using
# multiple impotance sampling
#  Modified by SPatel
#  Adapted by PTouloupou: trachoma model
### NOTES: using IUs instead of pixels
###
### Requires: Maps file, python code with parameter file, AMIS source
###############################################

get_scenario_id <- function(data, iscen) {
    Data = read.csv(data)
    scenar_group_pairs = unique(data.frame(Data$Scenario, Data$Group))

    ### TODO: Remove high prevalence group because...
    high_prevalence_group <- which(scenar_group_pairs$Data.Group == 7)
    scenar_group_pairs = scenar_group_pairs[-high_prevalence_group, ]

    return(scenar_group_pairs$Data.Scenario[iscen])
}

get_group_id <- function(data, iscen) {
    Data = read.csv(data)
    scenar_group_pairs = unique(data.frame(Data$Scenario, Data$Group))

    ### TODO: Remove high prevalence group because...
    high_prevalence_group <- which(scenar_group_pairs$Data.Group == 7)
    scenar_group_pairs = scenar_group_pairs[-high_prevalence_group, ]

    return(scenar_group_pairs$Data.Group[iscen])
}

sample_prevalence_map_at_IUs <- function(IU_indices, n.map.sampl, scenario_id) {
    prev = matrix(NA, ncol = n.map.sampl, nrow = length(IU_indices))
    sample_map <- function(IU_index) {
        set.seed(scenario_id) # For comparison with test data with `set.seed(Scen[iscen])`
        rnorm(n.map.sampl, Data$Logit[IU_index], sd = Data$Sds[IU_index])
    }
    L <- lapply(IU_indices, sample_map)
    prev <- sapply(L, function(x) exp(x)/(1+exp(x)))

    return(
        t(prev*100)
    )
}

iscen = 1

library(tmvtnorm)
library(mnormt)
library(mclust)

scenario_id <- get_scenario_id("./data/FinalDataPrev.csv", iscen)
group_id <- get_group_id("./data/FinalDataPrev.csv", iscen)

prefix <- sprintf("scen%g_group%g", scenario_id, group_id)
folder <- "output/"  # which folder to save final files to

Data = read.csv("./data/FinalDataPrev.csv")
IU_scen <- which(Data$Scenario == scenario_id & Data$Group == group_id)

prevalence_output <- sprintf("output/OutputPrev_scen%g_group%g.csv", scenario_id, group_id) # make sure this is consistent with main.py


source("AMIS_source.R")  # source code for AMIS

rprop0<-function(n){
  return(list(runif(n, min=0.05, max=0.175), runif(n, min=0, max=1)))
}

############## AMIS and MAP parameters ############
n.pixels<-length(IU_scen)  # Number of pixels OR IUs
n.map.sampl<-3000 # Number of samples for the map
ESS.R<-250 # Desired effective sample
delta<-5 # delta value (width for the Radon-Nikodym derivative) %
n.param<-2

T<-5; # max number of iterations
NN<-100  # Number of parameter sets in each iteration
N<-rep(NN,T)  # This allows to have different number of parameters sampled each iteration. Here it's the same  # different number of iterations might break code
#N[1] <- 50

prev <- sample_prevalence_map_at_IUs(IU_scen, n.map.sampl, scenario_id)
mean.prev<-sapply(1:n.pixels, function(a) mean(prev[a,]))

###################################################################
#          AMIS setup
####################################################################
# Set distribution for proposal: Student's t distribution
proposal=mvtComp(df=3); mixture=mclustMix();
dprop <- proposal$d
rprop <- proposal$r

param<-matrix(NA, ncol=n.param+1+1, nrow=sum(N))  # Matrix for parameter values, + prevalence and weights

IO_files_id <- function(t) {
    sprintf("scen%g_group%g_it%g", scenario_id, group_id, t)
}

###################################################################
#          Iteration 1.
####################################################################
t<-1  # Iteration
tmp<-rprop0(N[t])    #N[t] random draws of parameters from prior
x <- tmp[[1]]  # bet
y <- tmp[[2]]  # constant

seeds <- 1:N[t]
ans <- trachomAMIS::run_transmission_model(seeds, x, IO_files_id(1))

param[1:N[1],1]<-x
param[1:N[1],2]<-y
param[1:N[1],3]<-ans

first_weight <- rep(1, N[1])
WW <- trachomAMIS::compute_weight_matrix(prev, ans, delta, first_weight)
ess <- trachomAMIS::calculate_ess(WW)

cat( min(ess),  "", max(ess), "\n")

ESS<-matrix(ess, nrow=1, ncol=n.pixels)
### Copy variables for testing
### See tests/test_AMIS_trachoma.R
list_of_ESS <- list(ESS)
list_of_params <- list(param[1:N[1],1:3])

GG<-c()
Sigma <- list(NA, 10*T)
Mean<-list(NA, 10*T)
PP<-list(NA,T)
components <- list(GG = GG,
                   Sigma = Sigma,
                   Mean = Mean,
                   PP = PP)

set.seed(iscen)
for (t in 2:T) {

    WW <- update_according_to_ess_value(WW, ess, ESS.R)
    parameters <- param[1:sum(N[1:(t-1)]),1:2]
    clustMix <- trachomAMIS::evaluate_mixture(parameters, NN, WW, mixture)
    sampled_params <- trachomAMIS::sample_new_parameters(clustMix, N[t], rprop)
    param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),1]<-sampled_params$beta
    param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),2]<-sampled_params$constant

    seeds <- c((max(seeds)+1): (max(seeds)+N[t]))
    ans <-trachomAMIS::run_transmission_model(seeds, sampled_params$beta, IO_files_id(t))
    param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),3]<-ans
    components <- trachomAMIS::update_mixture_components(clustMix, components, t)
    first_weight <- trachomAMIS::compute_prior_proposal_ratio(components, t, T, N, beta = param[,1], constant = param[,2], dprop)

    all_sim_prevs<-param[1:sum(N[1:(t)]),3]
    WW <- trachomAMIS::compute_weight_matrix(prev, all_sim_prevs, delta, first_weight)
    ess <- trachomAMIS::calculate_ess(WW)

    cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))

    ESS<-rbind(ESS, as.numeric(ess))

    list_of_ESS[[t]] <- ESS
    print(sprintf("T = %g", t))
    param[1:5,]
    list_of_params[[t]] <- param[1:sum(N[1:t]),1:3]

    if(min(ess) >= ESS.R) break
}
