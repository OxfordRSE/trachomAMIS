devtools::load_all()

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

iscen <- 1
scenario_id <- 36
group_id <- 2
Data = read.csv("./data/FinalDataPrev.csv")
IU_scen <- which(Data$Scenario == scenario_id & Data$Group == group_id)
prev <- sample_prevalence_map_at_IUs(IU_scen, n.map.sampl = 3000, scenario_id)

############## Run AMIS ############
T <- 5
N<-rep(100,T)
param_and_weights <- trachomAMIS::amis(prevalence_map = prev, transmission_model = NULL, n_params = 2, nsamples = 100, IO_file_id = sprintf("scen%g_group%g",  scenario_id,  group_id), delta = 5, T = T, target_ess = 250)

source("tests/qa_ecdf_funcs.R")
ecdf_errors <- qa_plots(param_and_weights, prev, "tests/plots", Data$IUCodes[IU_scen], scenario_id, group_id)

for (irow in nrow(ecdf_errors)) {
    testthat::expect_lt(ecdf_errors$ecdf_max_distance[irow], 0.2)
    testthat::expect_lt(ecdf_errors$ecdf_sqrd_distance[irow], 0.05)
}
