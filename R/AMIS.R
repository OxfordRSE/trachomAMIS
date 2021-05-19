#' @export
amis <- function(prevalence_map, transmission_model, n_params, nsamples,
                 IO_file_id, delta = 5, T = 100, target_ess = 250) {
  param <- get_initial_parameters(nsamples)
  simulated_prevalences <- run_transmission_model(
    transmission_model,  
    seeds = 1:nsamples,
    parameters = param[, 1],
    IO_file_id
    )
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
  set.seed(iscen)
  seeds <- function(t) ((t - 1) * nsamples + 1):(t * nsamples)
  for (t in 2:T) {
    WW <- update_according_to_ess_value(WW, ess, target_ess)
    clustMix <- evaluate_mixture(param, nsamples, WW, mixture)
    new_params <- sample_new_parameters(clustMix, nsamples, prop$r)
    simulated_prevalences <- append(
      simulated_prevalences,
      run_transmission_model(transmission_model, seeds(t), new_params[, 1], IO_file_id)
    )
    components <- update_mixture_components(clustMix, components, t)
    param <- rbind(param, new_params)
    first_weight <- compute_prior_proposal_ratio(components, param, prop$d)
    WW <- compute_weight_matrix(prevalence_map, simulated_prevalences, delta, first_weight)
    ess <- calculate_ess(WW)
    if (min(ess) >= target_ess) break
  }

  allseeds <- 1:(T * nsamples)
  ret <- data.frame(allseeds, param[,-2], simulated_prevalences, t(WW))
  colnames(ret) <- c(
      "seeds",
      "beta",
      "sim_prev",
      sapply(1:dim(WW)[1], function(idx) sprintf("iu%g", idx))
      )
  return(ret)
}
