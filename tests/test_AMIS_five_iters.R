library(testthat)
devtools::load_all()

source("./AMIS_five_iterations.R")

### Test input file generated for the transmission model
### for first two iterations of AMIS

read_file <- function(iter, prefix) {
    filename <- paste(prefix, sprintf("_it%g.csv", iter), sep = "")
    read.csv(filename)
}

betas <- lapply(1:5, read_file, "./files/InputBet_scen36_group2")
expected_betas <- lapply(1:5, read_file, "./tests/test_data/InputBet_scen36_group2")

expect_equal(beta_iter_1, expected_beta_iter_1)
expect_equal(beta_iter_2, expected_beta_iter_2)

### Test Effective Sample Size (ESS) and content of `param` matrix
### for first two iterations of AMIS

read_matrix <- function(iter, prefix) {
    filename <- paste(prefix, sprintf("_iteration_%g.csv", iter), sep = "")
    matrix(scan(filename),
           nrow = dim(list_of_ESS[[iter]])[1],
           ncol = dim(list_of_ESS[[iter]])[2],
           byrow = TRUE)
}

expected_ESS <- lapply(1:5, read_matrix, "./tests/test_data/ESS")
expected_param <- lapply(1:5, read_matrix, "./tests/test_data/param")

t <- 1
expect_equal(expected_ESS[[1]], list_of_ESS[[t]])
expect_equal(expected_param[[1]], list_of_param[[t]][,1:3])

for (t in 1:5) {
    TRUE
}


