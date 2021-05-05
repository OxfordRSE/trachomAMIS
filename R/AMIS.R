#' @export
amis <- function(prevalence_map, transmission_model, n_params, nsamples, IO_file_id, delta = 5, T = 100, target_ess = 250) {
    param <- get_initial_parameters(nsamples)
    id <- function(t) paste(IO_file_id, sprintf("_it%g", t), sep = "")
    simulated_prevalences <- trachomAMIS::run_transmission_model(seeds = 1:nsamples, parameters = param[,1], id(1))
    WW <- trachomAMIS::compute_weight_matrix(prev, simulated_prevalences, delta, first_weight = rep(1, nsamples))
    ess <- trachomAMIS::calculate_ess(WW)
    cat( min(ess),  "", max(ess), "\n")

    components <- list(GG = c(0),
                       Sigma = list(),
                       Mean = list(),
                       PP = list())
                                        # Set distribution for proposal: Student's t distribution
    prop=mvtComp(df=3); mixture=mclustMix();
    set.seed(iscen)
    seeds <- function(t) c((t-1)*nsamples+1:t*nsamples)
    for (t in 2:T) {
        WW <- update_according_to_ess_value(WW, ess, target_ess)

        clustMix <- evaluate_mixture(param, nsamples, WW, mixture)
        param <- rbind(param, sample_new_parameters(clustMix, nsamples, prop$r))
        simulated_prevalences <- append(simulated_prevalences,
                                        run_transmission_model(seeds(t) , param[,1], id(t)))
        components <- trachomAMIS::update_mixture_components(clustMix, components, t)
        first_weight <- trachomAMIS::compute_prior_proposal_ratio(components, param, prop$d)
        WW <- trachomAMIS::compute_weight_matrix(prev, simulated_prevalences, delta, first_weight)
        ess <- trachomAMIS::calculate_ess(WW)
        cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))
        if(min(ess) >= target_ess) break
    }

    return(cbind(param, simulated_prevalences, deparse.level = 0))
}
