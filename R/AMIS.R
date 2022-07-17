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
  simulated_prevalences <- transmission_model(seeds = 1:nsamples, param[,1])
  weight_matrix <- compute_weight_matrix(
    prevalence_map,
    simulated_prevalences,
    amis_params,
    first_weight = rep(1-amis_params[["log"]], nsamples)
  )
  ess <- calculate_ess(weight_matrix,amis_params[["log"]])
  components <- list(
    G = c(0), # number of mixture components from proposal for each iteration (zero is for prior)
    Sigma = list(), # list of covariance matrices for each component
    Mean = list(), # list of means for each component
    probs = list() # probability of each component (unnormalised)
  )
  seeds <- function(t) ((t - 1) * nsamples + 1):(t * nsamples)
  niter <- 1
  for (t in 2:amis_params[["T"]]) {
    print(sprintf("AMIS iteration %g", t))
    mean_weights <- update_according_to_ess_value(weight_matrix, ess, amis_params[["target_ess"]],amis_params[["log"]]) 
    mixture <- weighted_mixture(param, amis_params[["mixture_samples"]], mean_weights, amis_params[["log"]])
    new_params <- sample_new_parameters(mixture, nsamples, amis_params[["df"]], prior, amis_params[["log"]])
    simulated_prevalences <- append(
      simulated_prevalences,
      transmission_model(seeds(t), new_params[,1])
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

  if(niter == amis_params[["T"]] && min(ess) < amis_params[["target_ess"]]) {
    msg <- sprintf(
      "Some locations did not reach target ESS (%g) after %g iterations",
      amis_params[["target_ess"]], niter
      )
    warning(msg)
  }

  allseeds <- 1:(niter * nsamples)
  ret <- data.frame(allseeds, param[,-2], simulated_prevalences, t(WW))
  colnames(ret) <- c(
      "seeds",
      "beta",
      "sim_prev",
      sapply(1:dim(WW)[1], function(idx) sprintf("iu%g", idx))
      )
  return(ret)
}
