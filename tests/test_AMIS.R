odevtools::load_all()

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

library(reticulate)

iscen <- 1
scenario_id <- get_scenario_id("./tests/test_data/FinalDataPrev.csv", iscen)
group_id <- get_group_id("./tests/test_data/FinalDataPrev.csv", iscen)

prev <- matrix(
    scan(file = "tests/test_data/prevalence_map.csv"),
    nrow = 3, byrow = T
)

reticulate::use_virtualenv("./.venv", required=TRUE)
model <- reticulate::import("trachoma")

############## Run AMIS ############
T <- 2
N<-rep(100,T)
param <- trachomAMIS::amis(prevalence_map = prev, transmission_model = model$Trachoma_Simulation, n_params = 2, nsamples = 100, IO_file_id = sprintf("scen%g_group%g",  scenario_id,  group_id), delta = 5, T = T, target_ess = 250, mda_file = "files/InputMDA_scen36_group2.csv", jobid = 1, seed = 1)

expected_param <- read.csv("./tests/test_data/param_iteration_2.csv")
testthat::expect_equal(param, expected_param)
