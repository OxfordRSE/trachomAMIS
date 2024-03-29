#' Run the AMIS algorithm to fit a transmission model to a map
#'
#' For details of the algorithm, see
#' Integrating geostatistical maps and infectious disease transmission models
#' using adaptive multiple importance sampling.
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' Ann. Appl. Stat. 15 (4) 1980 - 1998, December 2021.
#' DOI: https://doi.org/10.1214/21-AOAS1486
#'
#' @param prevalence_map An L x M matrix containing samples from the fitted prevalence map, where L is the number of locations and M the number of samples.
#' The location names are inherited from \code{rownames(prevalence_map)} if possible. Alternatively, a list with one entry for each timepoint.
#' Each entry must be a list containing objects \code{data} (an L x M matrix of data as above);
#' and \code{likelihood} a function taking arguments \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.
#' If a likelihood is not specified then it is assumed that
#' the data consist of samples from a geo-statistical model and empirical methods are used.
#' @param transmission_model A function taking a vector of n seeds and an n x d matrix of parameter vectors as inputs
#'  and producing a n x timepoints MATRIX of prevalences as output (it must be a matrix even when timepoints==1).
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
#' \item{\code{log} logical indicating if calculations are to be performed on log scale.}
#' \item{\code{max_iters}maximum number of AMIS iterations.}
#' \item{\code{breaks}optional vector specifying the breaks for the histogram.}
#' }
#' @param seed Optional seed for the random number generator
#' @return A dataframe of the sampled parameters, simulation seed, and weight in each location.
#' @export
amis <- function(prevalence_map, transmission_model, prior, amis_params, seed = NULL) {
  if (is.matrix(prevalence_map) || is.data.frame(prevalence_map)) {prevalence_map=list(list(data=prevalence_map))}
  # add some checks to beginning of function with helpful error messages?
  if(!is.null(seed)) set.seed(seed)
  nsamples <- amis_params[["nsamples"]]
  print("AMIS iteration 1")
  # Sample first set of parameters from the prior
  param <- prior$rprior(nsamples)
  # to avoid duplication, evaluate prior density now.
  prior_density<-sapply(1:nsamples,function(b) {prior$dprior(param[b,],log=amis_params[["log"]])})
  # Simulate from transmission model
  simulated_prevalences <- transmission_model(seeds = 1:nsamples, param)
  # to avoid duplication, evaluate likelihood now.
  likelihoods<- compute_likelihood(prevalence_map,simulated_prevalences,amis_params)
  weight_matrix <- compute_weight_matrix(
    likelihoods,
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
    if ((amis_params[["log"]] && max(mean_weights)==-Inf) || (!amis_params[["log"]] && max(mean_weights)==0)) {stop("No weight on any particles for locations in the active set.\n")}
    mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
    components <- update_mixture_components(mixture, components, t)
    new_params <- sample_new_parameters(mixture, nsamples, amis_params[["df"]], prior, amis_params[["log"]])
    param <- rbind(param, new_params$params)
    prior_density <- c(prior_density,new_params$prior_density)
    new_prevalences <- transmission_model(seeds(t), new_params$params)
    simulated_prevalences <- rbind(simulated_prevalences,new_prevalences)
    likelihoods <- compute_likelihood(prevalence_map,new_prevalences,amis_params,likelihoods)
    first_weight <- compute_prior_proposal_ratio(components, param, prior_density, amis_params[["df"]], amis_params[["log"]]) # Prior/proposal
    weight_matrix <- compute_weight_matrix(likelihoods, simulated_prevalences, amis_params, first_weight) # RN derivative (shd take all amis_params)
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
  ret <- data.frame(allseeds, param, simulated_prevalences, weight_matrix)
  if (is.null(rownames(prevalence_map[[1]]$data))) {
    iunames<-sapply(1:dim(weight_matrix)[2], function(idx) sprintf("iu%g", idx))
  } else {
    iunames<-rownames(prevalence_map[[1]]$data)
  }
  if (is.null(colnames(param))) {
    paramnames<-sapply(1:dim(param)[2], function(idx) sprintf("param%g", idx))
  } else {
    paramnames<-colnames(param)
  }
  if (is.null(colnames(simulated_prevalences))) {
    prevnames<-sapply(1:dim(simulated_prevalences)[2], function(idx) sprintf("prev%g", idx))
  } else {
    prevnames<-paste0("prev",colnames(simulated_prevalences))
  }
  colnames(ret) <- c("seeds",paramnames,prevnames,iunames)
  return(ret)
}
