devtools::load_all()

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

library(reticulate)

scenario_id <- get_scenario_id("./data/FinalDataPrev.csv", iscen)
group_id <- get_group_id("./data/FinalDataPrev.csv", iscen)

prefix <- sprintf("scen%g_group%g", scenario_id, group_id)
folder <- "output/"  # which folder to save final files to

Data = read.csv("./data/FinalDataPrev.csv")
IU_scen <- which(Data$Scenario == scenario_id & Data$Group == group_id)

prevalence_output <- sprintf("output/OutputPrev_scen%g_group%g.csv", scenario_id, group_id) # make sure this is consistent with main.py

prev <- sample_prevalence_map_at_IUs(IU_scen, n.map.sampl = 3000, scenario_id)

reticulate::use_virtualenv("./.venv", required=TRUE)
model <- reticulate::import("trachoma")

############## Run AMIS ############
T <- 5
N<-rep(100,T)
param <- trachomAMIS::amis(prevalence_map = prev, transmission_model = model$Trachoma_Simulation, n_params = 2, nsamples = 100, IO_file_id = sprintf("scen%g_group%g",  scenario_id,  group_id), delta = 5, T = T, target_ess = 250)

expected_param <- matrix(scan("./tests/test_data/param_iteration_5.csv"),
                         nrow = 500,
                         ncol = 3,
                         byrow = TRUE)
testthat::expect_equal(param, expected_param)
