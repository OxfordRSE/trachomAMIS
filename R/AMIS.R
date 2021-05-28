#' @export
amis <- function(prevalence_map, transmission_model, nsamples,
                 IO_file_id, delta = 5, T = 100, target_ess = 250,
                 mda_file, seed = NULL) {

  if(!is.null(seed)) set.seed(seed)
    
  print("AMIS iteration 1")
  param <- get_initial_parameters(nsamples)
  simulated_prevalences <- transmission_model(seeds = 1:nsamples, param[,1])
  WW <- compute_weight_matrix(
    prevalence_map,
    simulated_prevalences,
    delta,
    first_weight = rep(1, nsamples)
  )
  ess <- calculate_ess(WW)
  components <- list(
    GG = c(0),
    Sigma = list(),
    Mean = list(),
    PP = list()
  )
  prop <- mvtComp(df = 3)
  mixture <- mclustMix()
  seeds <- function(t) ((t - 1) * nsamples + 1):(t * nsamples)
  niter <- 1
  for (t in 2:T) {
    print(sprintf("AMIS iteration %g", t))
    WW <- update_according_to_ess_value(WW, ess, target_ess)
    clustMix <- evaluate_mixture(param, nsamples, WW, mixture)
    new_params <- sample_new_parameters(clustMix, nsamples, prop$r)
    simulated_prevalences <- append(
      simulated_prevalences,
      transmission_model(seeds(t), new_params[,1])
    )
    components <- update_mixture_components(clustMix, components, t)
    param <- rbind(param, new_params)
    first_weight <- compute_prior_proposal_ratio(components, param, prop$d)
    WW <- compute_weight_matrix(prevalence_map, simulated_prevalences, delta, first_weight)
    ess <- calculate_ess(WW)
    niter <- niter + 1
    if (min(ess) >= target_ess) break
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
