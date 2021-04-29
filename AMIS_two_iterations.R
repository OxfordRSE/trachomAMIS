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

# function to draw from prior
dprop0<-function(a,b){
  return(dunif(a, min=0.05, max=0.175)*dunif(b, min=0, max=1))
}
rprop0<-function(n){
  return(list(runif(n, min=0.05, max=0.175), runif(n, min=0, max=1)))
}

############## AMIS and MAP parameters ############
n.pixels<-length(IU_scen)  # Number of pixels OR IUs
n.map.sampl<-3000 # Number of samples for the map
ESS.R<-250 # Desired effective sample
delta<-5 # delta value (width for the Radon-Nikodym derivative) %
n.param<-2

T<-100; # max number of iterations
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
Sigma <- list(NA, 10*T)
Mean<-list(NA, 10*T)
PP<-list(NA,T)
GG<-list(NA,T)


###################################################################
#          Iteration 1.
####################################################################
t<-1  # Iteration
tmp<-rprop0(N[t])    #N[t] random draws of parameters from prior
x <- tmp[[1]]  # bet
y <- tmp[[2]]  # constant
seed <- c(1:N[t])
allseed <- seed
input_params <- cbind(seed, x)
colnames(input_params) = c("randomgen", "bet")

inputbeta <- sprintf("files/InputBet_scen%g_group%g_it1.csv", scenario_id, group_id)
write.csv(input_params, file=inputbeta, row.names=FALSE)

### Run Python

### read in python output file
output_file <- sprintf("output/OutputPrev_scen%g_group%g_it1.csv", scenario_id, group_id)
ans <- trachomAMIS::read_simulated_prevalence(output_file)

w<-sapply(1:length(ans), function(i) length(which((prev>ans[i]-delta/2) &(prev<=ans[i]+delta/2)))/length(which((ans>ans[i]-delta/2) & (ans<=ans[i]+delta/2))))   #weights over all IUs

param[1:N[1],1]<-x
param[1:N[1],2]<-y
param[1:N[1],3]<-ans
param[1:N[1],4]<- w

prop<-param[1:N[1],]

first_weight <- rep(1, N[1])
WW <- trachomAMIS::compute_weight_matrix(prev, ans, delta, first_weight)
ess <- trachomAMIS::calculate_ess(WW)

cat( min(ess),  "", max(ess), "\n")

ESS<-matrix(ess, nrow=1, ncol=n.pixels)
### Copy variables for testing
### See tests/test_AMIS_trachoma.R
ESS_iteration_1 <- ESS
param_iteration_1 <- param[1:N[1],]



###################################################################
#          Iteration 2+
####################################################################

set.seed(iscen)
stop<-0
  
t<-t+1
cat(c("Iteration: ", t,", min(ESS): ", min(ess),"\n"))

wh<-which(ess>=ESS.R)
W1<-WW; W1[wh,]<-0

w1<- c(colSums(W1))


J<-sample(1:sum(N[1:(t-1)]), NN, prob= w1, replace=T)
xx<-param[J,1:2]
clustMix <- mixture(xx)

G <- clustMix$G
cluster <- clustMix$cluster
sampled_params <- trachomAMIS::sample_new_parameters(clustMix, N[t])

### Components of the mixture
ppt <- clustMix$alpha
muHatt <- clustMix$muHat
varHatt <- clustMix$SigmaHat
GG[[t-1]]<-G
G1<-0; G2<-G
if(t>2) {
    G1<-sum(sapply(1:(t-2), function(a) GG[[a]]))
    G2<-sum(sapply(1:(t-1), function(a) GG[[a]]))
}
for(i in 1:G){
    Sigma[[i+G1]] <- varHatt[,,i]
    Mean[[i+G1]] <- muHatt[i,]
    PP[[i+G1]]<-ppt[i]   ### scale by number of points
}


print("done sampling")
print(Sys.time())

seed <- c((max(seed)+1): (max(seed)+N[t]))
allseed <- c(allseed, seed)
input_params <- cbind(seed, sampled_params$beta)
colnames(input_params) = c("randomgen", "bet")

inputbeta <- sprintf("files/InputBet_scen%g_group%g_it2.csv", scenario_id, group_id)
write.csv(input_params, file=inputbeta, row.names=FALSE)

prevalence_output <- sprintf("output/OutputPrev_scen%g_group%g_it2.csv", scenario_id, group_id)
res <- read.csv(prevalence_output) # read python output file
ans <- 100*res[,dim(res)[2]]

param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),1]<-sampled_params$beta
param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),2]<-sampled_params$constant
param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),3]<-ans

prop.val <- sapply(1:sum(N[1:t]),function(b)  sum(sapply(1:G2, function(g) PP[[g]] * dprop(param[b,1:2],mu= Mean[[g]], Sig=Sigma[[g]]))) + dprop0(param[b,1], param[b,2]))   ## FIX to be just the proposal density ALSO scale by number of points

first_weight <- sapply(1:sum(N[1:t]), function(b) dprop0(param[b,1], param[b,2])/prop.val[b])   # prior/proposal


ans<-param[1:sum(N[1:(t)]),3]

ess<-c()
WW<-matrix(NA, nrow=n.pixels, ncol=sum(N[1:(t)]))
for(i in 1:n.pixels){
    w<-sapply(1:length(ans), function(j) length(which((prev[i,]>ans[j]-delta/2) &(prev[i,]<=ans[j]+delta/2)))/sum(first_weight[which((ans>ans[j]-delta/2) & (ans<=ans[j]+delta/2))]) )   # f/g from AMIS paper #### FIX
    ww<-w*first_weight; ww<-ww/sum(ww)  # second weighting, normalizing
    WW[i,]<-ww
    www<-(sum((ww)^2))^(-1)
    ess<-c(ess, www)
}

cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))

ESS<-rbind(ESS, as.numeric(ess))

w1<-c(colSums(WW))
param[1:sum(N[1:(t)]),4]<-w1

### Copy variables for testing
### See tests/test_AMIS_trachoma.R
ESS_iteration_2 <- ESS
param_iteration_2 <- param[1:N[1],]
