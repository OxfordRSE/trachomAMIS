dprop0 <- function(x) {
  return(dunif(x[1], min = 0.05, max = 0.175) * dunif(x[2], min = 0, max = 1))
}

#' Sample initial parameters from prior
#'
#' Sample N initial values for two parameters, uniformly distributed
#' between MIN and MAX.
#'
#' @param n The number of samples to draw for each parameter.
#' @param min The lower bound for the uniform distribution.
#' @param max The upper bound for the uniform distribution.
#' @return A Nx2 matrix
#' @examples
#' get_initial_parameters(100, 0.05, 0.175)
get_initial_parameters <- function(n) {
  ## Samples initial parameter from uniform distrib (prior)
  init_beta_samples <- runif(n, min = 0.05, max = 0.175)
  init_constant_samples <- runif(n, min = 0, max = 1)
  return(
    matrix(c(init_beta_samples, init_constant_samples),
      ncol = 2
    )
  )
}

#' Write input file required by transmission model
#'
#' Write a 2 columns csv file named INPUT_FILE to disk. The first column SEEDS
#' contains the seed for each parameter value. Second column BETA contains
#' corresponding parameter values.
#'
#' @param seeds A vector containing seed value for each parameter sample
#'     (double)
#' @param beta A vector containing samples of the beta parameter (double)
#' @param input_file The name of the input file (string)
write_model_input <- function(seeds, beta, input_file) {
  input_params <- cbind(seeds, beta)
  colnames(input_params) <- c("randomgen", "bet")
  write.csv(input_params, file = input_file, row.names = FALSE)
}

run_transmission_model <- function(model_func, seeds, parameters, id) {
  input_file <- paste("files/InputBet_", id, ".csv", sep = "")
  write_model_input(seeds, parameters, input_file)
  output_file <- paste("output/OutputPrev_", id, ".csv", sep = "")
  inputMDA <- paste("files/InputMDA_", id, ".csv", sep = "")
  infect_output <- paste("output/InfectFilePath_", id, ".csv", sep = "")
  model_func(input_file, inputMDA, output_file, infect_output,
  model_func(input_file, mda_file, output_file, infect_output,
             SaveOutput = F, OutSimFilePath = NULL, InSimFilePath = NULL)
  res <- read.csv(output_file)
  return(100 * res[, dim(res)[2]])
}

#' Compute weight matrix
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' Implementation Unit (IU). One row per sample, one column per IU.  Each weight
#' is computed based on the Randon-Nikodim derivative, taking into account
#' geostatistical prevalence data for the specific IU and the prevalence value
#' computed from the transmission model for the specific parameter sample.
#'
#' @param prev_data The geostatistical prevalence data. A n x m matrix where n
#'     is the number of IUs and m the number of prevalence samples. (double)
#' @param prev_sim A vector containing the siumulated prevalence value for each
#'     parameter sample. (double)
#' @param delta A number. (double)
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as \code{prev_sim}.
compute_weight_matrix <- function(prev_data, prev_sim, delta, first_weight) {
  n_IUs <- dim(prev_data)[1]
  weight_mat <- matrix(NA, nrow = n_IUs, ncol = length(prev_sim))

  radon_niko_deriv <- function(idx, prev_data_for_IU) {
    f <- length(which((prev_data_for_IU > prev_sim[idx] - delta / 2) & (prev_data_for_IU <= prev_sim[idx] + delta / 2)))
    g <- sum(first_weight[which((prev_sim > prev_sim[idx] - delta / 2) & (prev_sim <= prev_sim[idx] + delta / 2))])

    return(f / g)
  }

  for (i in 1:n_IUs) {
    w <- sapply(1:length(prev_sim), radon_niko_deriv, prev_data[i, ])
    w <- w * first_weight
    if (sum(w) > 0) w <- w / sum(w)
    weight_mat[i, ] <- w
  }
  return(weight_mat)
}

#' Compute the current effective sample size
#'
#' This function returns the effective sample size (ESS) for each Implementation
#' Unit (IU). For each column of the weight matrix WEIGHT_MAT, the ESS component
#' is computed as
#'
#' \deqn{\left(\sum_{i=1}^{N} w_{i}^2 \right)^{-1}},
#'
#' where \eqn{$N$} is the number of sampled parameter values for each IU.
#'
#' @param weight_mat The weight matrix. A n x m matrix with n the number of IUS
#'     and m the number of sampled parameter values.
#' @return A vector containing the ESS value for each IU.
#'
#' @seealso \code{\link{compute_weight_matrix}}
calculate_ess <- function(weight_mat) {
  ess_for_IU <- function(weights_for_IU) {
    if (sum(weights_for_IU) == 0) {
      return(0)
    }
    return(
      (sum((weights_for_IU)^2))^(-1)
    )
  }

  return(
    apply(weight_mat, 1, ess_for_IU)
  )
}

#' Update weight matrix
#'
#' This function sets to 0 the rows of the weight matrix WEIGHT_MATRIX for which
#' the effective sample size ESS is above a target size TARGET_SIZE.
#'
#' @param weight_matrix The weight_matrix as returned by
#'     \link{compute_weight_matrix}
#' @param ess The effective sample size vector as returned by
#'     \link{calculate_ess}
#' @param target_size A number representing the target size for the sample.
#' @return The updated weight matrix. Its size is unchanged.
update_according_to_ess_value <- function(weight_matrix, ess, target_size) {
  rows_to_nullify <- which(ess >= target_size)
  weight_matrix[rows_to_nullify, ] <- 0
  return(weight_matrix)
}

#' Evaluate the mixture function as returned by \link{mclustMix}
#'
#' Evaluate the mixture function MIXTURE for an subset of NSAMPLES parameters
#' values sampled according to the weight matrix WEIGHT_MATRIX.
#' @param parameters A matrix m x 2 matrix containing the sampled values for the
#'     2 parameters.
#' @param nsamples The size of the parameter subset to evaluate the mixture
#'     function from (integer)
#' @param weight_matrix The weight_matrix as returned by \link{compute_weight_matrix}
#' @param mixture The mixture function as returned by function \link{mclustMix}
#' @return A list of the mixture components (see function \code{\link{mclustMix}})
#'     \describe{
#'       \item{\code{alpha}}{The probability?}
#'       \item{\code{muHat}}{The mean?}
#'       \item{\code{SigmaHat}}{The variance?}
#'       \item{\code{G}}{??}
#'       \item{\code{cluster}}{??}
#'     }
#'
#' @seealso \code{\link{mclustMix}}
evaluate_mixture <- function(parameters, nsamples, weight_matrix, mixture) {
  sampled_idx <- sample(
    1:dim(parameters)[1],
    nsamples,
    prob = colSums(weight_matrix),
    replace = T
  )
  return(
    mixture(parameters[sampled_idx, ])
  )
}

#' Sample new parameters
#'
#' This function generates NSAMPLES new model parameter values according the
#' proposal distribution RPROP and the mixture components CLUSTMIX.
#'
#' @param clustMix A list of mixture components as returned by
#'     \code{\link{evaluate_mixture}}
#' @param nsamples A number of new parameters to sample (integer)
#' @param rprop The proposal distribution as returned by \code{mvtComp()$r}.
#' @return A NSAMPLES x 2 matrix containing the sampled parameter values.
#'
#' @seealso \code{\link{mvtComp}}
sample_new_parameters <- function(clustMix, n_samples, rprop) {
  x <- c()
  y <- c()
  while (length(x) < n_samples) {
    compo <- sample(1:clustMix$G, 1, prob = clustMix$alpha)
    x1 <- t(
      rprop(1, clustMix$muHat[compo, ], clustMix$SigmaHat[, , compo])
    )
    new.param <- as.numeric(x1)
    if (dprop0(new.param) > 0) {
      x <- c(x, new.param[1])
      y <- c(y, new.param[2])
    }
  }
  return(
    matrix(c(x, y), ncol = 2)
  )
}

#' Update the components of the mixture
#'
#' This function updates the components of the mixture COMPONENTS according to
#' the current mixture CLUSTMIX generated at iteration T.
#'
#' @param clustMix A list of mixture components as returned by
#'     \code{\link{evaluate_mixure}}
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{GG}}{A numeric vector}
#'       \item{\code{Sigma}}{A list}
#'       \item{\code{Mean}}{A list}
#'       \item{\code{PP}}{A list}
#'     }
#' @param t The current iteration index (integer)
#' @return The updated \code{components} list
#'
#' @seealso \code{\link{evaluate_mixture}}, \code{\link{mclustMix}}
update_mixture_components <- function(clustMix, components, t) {
  components$GG[t] <- clustMix$G
  G1 <- sum(components$GG[1:(t - 1)])
  for (i in 1:clustMix$G) {
    components$Sigma[[i + G1]] <- clustMix$SigmaHat[, , i]
    components$Mean[[i + G1]] <- clustMix$muHat[i, ]
    components$PP[[i + G1]] <- clustMix$alpha[i] ### scale by number of points
  }
  return(components)
}


#' Compute the prior/proposal ratio
#'
#' This function returns the ratio between the prior and proposal distribution
#' for each sampled parameter value (i.e. each row in PARAM).
#'
#' This function returns
#' \deqn(\frac{\pi(\theta _j)}{\frac{1}{N_1 + ... + N_t}\sum_{u=1}^{t}N_u \phi_u(\theta _j)})
#' See step (4) if the AMIS algorithm in
#'
#' Integrating Geostatistical Maps And Transmission Models Using Adaptive
#' Multiple Importance Sampling
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' medRxiv 2020.08.03.20146241;
#' doi: \url{https://doi.org/10.1101/2020.08.03.20146241}
#'
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{GG}}{A numeric vector}
#'       \item{\code{Sigma}}{A list}
#'       \item{\code{Mean}}{A list}
#'       \item{\code{PP}}{A list}
#'     }
#' @param param A 2 column matrix containing the values for the two parameters.
#' @param dprop The proposal distribution as returned by \code{mvtComp()$d}.
#' @return A vector contaning the prior/proposal ratio for each row in
#'     \code{param}
compute_prior_proposal_ratio <- function(components, param, dprop) {
  PP <- components$PP
  Sigma <- components$Sigma
  Mean <- components$Mean

  G2 <- sum(components$GG)
  prop.val <- sapply(1:nrow(param), function(b) sum(sapply(1:G2, function(g) PP[[g]] * dprop(param[b, ], mu = Mean[[g]], Sig = Sigma[[g]]))) + dprop0(param[b, ])) ## FIX to be just the proposal density ALSO scale by number of points

  first_weight <- sapply(1:nrow(param), function(b) dprop0(param[b, ]) / prop.val[b]) # prior/proposal

  return(first_weight)
}
