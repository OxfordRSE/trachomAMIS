library(testthat)

source("./AMIS_two_iterations.R")

### Test input file generated for the transmission model
### for first two iterations of AMIS

beta_iter_1 <- read.csv("./files/InputBet_scen36_group2_it1.csv")
beta_iter_2 <- read.csv("./files/InputBet_scen36_group2_it2.csv")
expected_beta_iter_1 <- read.csv("./tests/test_data/InputBet_scen36_group2_it1.csv")
expected_beta_iter_2 <- read.csv("./tests/test_data/InputBet_scen36_group2_it2.csv")

expect_equal(beta_iter_1, expected_beta_iter_1)
expect_equal(beta_iter_2, expected_beta_iter_2)

### Test Effective Sample Size (ESS) and content of `param` matrix
### for first two iterations of AMIS

expected_ESS_iter_1 <- matrix(scan("./tests/test_data/ESS_iteration_1.csv"),
                              nrow = dim(ESS_iteration_1)[1],
                              ncol = dim(ESS_iteration_1)[2],
                              byrow = TRUE)
expected_ESS_iter_2 <- matrix(scan("./tests/test_data/ESS_iteration_2.csv"),
                              nrow = dim(ESS_iteration_2)[1],
                              ncol = dim(ESS_iteration_2)[2],
                              byrow = TRUE)

expected_param_iter_1 <- matrix(scan("./tests/test_data/param_iteration_1.csv"),
                              nrow = dim(param_iteration_1)[1],
                              ncol = dim(param_iteration_1)[2],
                              byrow = TRUE)
expected_param_iter_2 <- matrix(scan("./tests/test_data/param_iteration_2.csv"),
                              nrow = dim(param_iteration_2)[1],
                              ncol = dim(param_iteration_2)[2],
                              byrow = TRUE)

expect_equal(ESS_iteration_1, expected_ESS_iter_1)
expect_equal(ESS_iteration_2, expected_ESS_iter_2)
expect_equal(param_iteration_1, expected_param_iter_1)
expect_equal(param_iteration_2, expected_param_iter_2)


