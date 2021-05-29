devtools::load_all()

args <- commandArgs(trailingOnly = T)
if (!(length(args) == 0)) {
    plot_dir <- args[1]
} else {
    plot_dir <- plot_dir <- file.path("tests", "ecdf_plots")
}
if (!dir.exists(plot_dir)) dir.create(plot_dir)

reticulate::use_virtualenv("./.venv", required=TRUE)
module <- reticulate::import("trachoma")
run_model <- module$Trachoma_Simulation

make_file_path <- function(prefix) {
    file.path("model_io", sprintf("%s_scen36.csv", prefix))
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

prev <- matrix(
    scan(file = "tests/test_data/prevalence_map.csv"),
    nrow = 3, byrow = T
)

############## Run AMIS ############
params <- list("nsamples" = 100,
               "T" = 100,
               "target_ess" = 250,
               "delta" = 5)
param <- trachomAMIS::amis(prev, wrapped_model, params)

## Codes for IUs in group 2 with scenario 36
iucodes <- c("ETH18551", "ETH18644", "ETH18541")

source("tests/qa_ecdf_funcs.R")
ecdf_errors <- qa_plots(param, prev, plot_dir, iucodes, scenid = 36, grpid = 2)
