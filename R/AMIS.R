# SS note: Provide option to supply prevalence_map as an analytical function (likelihood and data), rather than Monte Carlo samples
# SS note: Also implement Renata's RN derivative?
# SS note: Also implement minimum error RN derivative?

#' Run the AMIS algorithm to fit a transmission model to a map
#'
#' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map, where L is the number of locations and M the number of samples.
#' The location names are inherited from \code{rownames(prevalence_map)} if possible. If instead, an analytic form of the prevalence map is available,
#' supply a list containing objects \code{data} (an L x M matrix of data) and \code{likelihood} a function taking a row of data and the output from the transmission
#' model as arguments (and logical \code{log}) and returning the (log)-likelihood. 
#' @param transmission_model A function taking a vector of seeds and a matrix of parameter vectors as inputs.
#' @param prior A list containing the functions \code{dprior} and \code{rprior} (density and RNG).
#' \code{rprior} must produce an n by d MATRIX of parameters, even when d=1.
#' parameter names are inherited from the \code{colnames} from the output of \code{rprior} if possible.
#' @param amis_params A list containing the control parameters for the AMIS algorithm
#' \describe{
#' \item{\code{delta}the smoothing parameter in the empirical RN derivative (usually 0.01).}
#' \item{\code{nsamples}the number of new samples drawn within each AMIS iteration.}
#' \item{\code{mixture_samples}the number of samples used to represent the weighted parameters in the mixture fitting.}
#' \item{\code{df}the degrees of freedom in the t-distributions, used to yield a heavy tailed proposal.}
#' \item{\code{target_ess}the target effective sample size.}
#' \item{\code{log} logicalindicating if calculations are to be performed on log scale.} 
#' \item{\code{max_iters}maximum number of AMIS iterations.}
#' \item{\code{method}string describing the appropriate method to use, e.g. empirical.}
#' }
#' @param seed Optional seed for the random number generator
#' @return A list containing a dataframe of the sampled parameters, simulation seed, and weight in each location, plus a vector
#' #' called ess containing the obtained ess at each location.  
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params, seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  nsamples <- amis_params[["nsamples"]]
  print("AMIS iteration 1")
  # Sample first set of parameters from the prior
  param <- prior$rprior(nsamples)
  # to avoid duplication, evaluate prior density now.
  prior_density<-sapply(1:nsamples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
  # Simulate from transmission model
  simulated_prevalences <- transmission_model(seeds = 1:nsamples, param) # SS changed to pass ALL parameters to the transmission model
  weight_matrix <- compute_weight_matrix(
    prevalence_map,
    simulated_prevalences,
    amis_params,
    first_weight = rep(1-amis_params[["log"]], nsamples)
  )
  ess <- calculate_ess(weight_matrix,amis_params[["log"]])
  # Make object to store the components of the AMIS mixture.
  components <- list(
    G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
    Sigma = list(), # list of covariance matrices for each component
    Mean = list(), # list of means for each component
    probs = list() # probability of each component (unnormalised)
  )
  seeds <- function(t) ((t - 1) * nsamples + 1):(t * nsamples)  #function to calculate the seeds for iteration t.
  niter <- 1 # number of completed iterations
  for (t in 2:amis_params[["max_iters"]]) {
    print(sprintf("AMIS iteration %g", t))
    mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]]) 
    mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
    new_params <- sample_new_parameters(mixture, nsamples, amis_params[["df"]], prior, amis_params[["log"]])
    simulated_prevalences <- append(
      simulated_prevalences,
      transmission_model(seeds(t), new_params$params)
    )
    components <- update_mixture_components(mixture, components, t)
    param <- rbind(param, new_params$params)
    prior_density <- c(prior_density,new_params$prior_density)
    first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
    weight_matrix <- compute_weight_matrix(prevalence_map, simulated_prevalences, amis_params, first_weight) # RN derivative (shd take all amis_params)
    ess <- calculate_ess(weight_matrix,amis_params[["log"]])
    niter <- niter + 1
    if (min(ess) >= amis_params[["target_ess"]]) break
  }

  if(niter == amis_params[["max_iters"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg)
  }
  allseeds <- 1:(niter * nsamples)
  ret <- data.frame(allseeds, param, simulated_prevalences, t(weight_matrix))
  if (is.null(rownames(prevalence_map))) {
    iunames<-sapply(1:dim(weight_matrix)[1], function(idx) sprintf("iu%g", idx))
  } else {
    iunames<-rownames(prevalence_map)
  }
  if (is.null(colnames(param))) {
    paramnames<-sapply(1:dim(param)[2], function(idx) sprintf("param%g", idx))
  } else {
    paramnames<-colnames(param)
  }
  colnames(ret) <- c("seeds",paramnames,"sim_prev",iunames)
  return(ret)
}
