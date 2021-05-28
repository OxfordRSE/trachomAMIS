devtools::load_all()
library(reticulate)

scenario_id <- 36; group_id <- 2

prev <- matrix(
    scan(file = "tests/test_data/prevalence_map.csv"),
    nrow = 3, byrow = T
)

reticulate::use_virtualenv("./.venv", required=TRUE)
module <- reticulate::import("trachoma")
run_model <- module$Trachoma_Simulation

make_file_path <- function(prefix) {
    file.path("model_io", sprintf("%s_scen%g.csv", prefix, scenario_id))
}
input_file <- make_file_path("InputBet")
output_file <- make_file_path("OutputPrev")
infect_output <- make_file_path("InfectOutput")
mda_file <- "./tests/test_data/InputMDA_scen36_group2.csv"

wrapped_model <- function(seeds, parameters) {
  ## write input on disk
  dir.create("model_io")
  on.exit(unlink("model_io", recursive = TRUE))
  write.csv(cbind(seeds, parameters), input_file, row.names = F)
  ## run model
  run_model(input_file, mda_file, output_file, infect_output,
            SaveOutput = F, OutSimFilePath = NULL,
            InSimFilePath = NULL)
  ## read results and return obsvervable
  res <- read.csv(output_file)
  return(100 * res[, dim(res)[2]])
}

############## Run AMIS ############
params <- list("nsamples" = 100,
               "T" = 2,
               "target_ess" = 250,
               "delta" = 5)
param <- trachomAMIS::amis(prev, wrapped_model, params, seed = 1)

expected_param <- read.csv("./tests/test_data/param_iteration_2.csv")
testthat::expect_equal(param, expected_param)
