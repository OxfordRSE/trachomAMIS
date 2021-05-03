#' @export
amis <- function(prevalence_map, transmission_model, n_params, N, IO_file_id, delta = 5, T = 100, target_ess = 250) {
    param<-matrix(NA, ncol=n_param+1+1, nrow=sum(N))  # Matrix for parameter values, + prevalence and weights
    tmp<-rprop0(N[1])    #N[t] random draws of parameters from prior
    x <- tmp[[1]]  # bet
    y <- tmp[[2]]  # constant

    sim_prev <- trachomAMIS::run_transmission_model(seeds = 1:N[1], x, IO_files_id)
    param[1:N[1],1]<-x
    param[1:N[1],2]<-y
    param[1:N[1],3]<-ans

    WW <- trachomAMIS::compute_weight_matrix(prev, ans, delta, first_weight = rep(1, N[1]))
    ess <- trachomAMIS::calculate_ess(WW)
    cat( min(ess),  "", max(ess), "\n")

    GG<-list(NA,T)
    Sigma <- list(NA, 10*T)
    Mean<-list(NA, 10*T)
    PP<-list(NA,T)

    # Set distribution for proposal: Student's t distribution
    proposal=mvtComp(df=3); mixture=mclustMix();
    dprop <- proposal$d
    rprop <- proposal$r

    set.seed(iscen)
    for (t in 2:T) {
        WW <- update_according_to_ess_value(WW, ess, ESS.R)
        parameters <- param[1:sum(N[1:(t-1)]),1:2]
        clustMix <- trachomAMIS::evaluate_mixture(parameters, NN, WW, mixture)
        sampled_params <- trachomAMIS::sample_new_parameters(clustMix, N[t])
        param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),1]<-sampled_params$beta
        param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),2]<-sampled_params$constant
        seeds <- c((max(seeds)+1): (max(seeds)+N[t]))
        ans <-trachomAMIS::run_transmission_model(seeds, sampled_params$beta, IO_files_id)
        param[(sum(N[1:(t-1)])+1):sum(N[1:(t)]),3]<-ans
        first_weight <- trachomAMIS::compute_prior_proposal_ratio(clustMix, t, T, N, beta = param[,1], constant = param[,2])
        all_sim_prevs<-param[1:sum(N[1:(t)]),3]
        WW <- trachomAMIS::compute_weight_matrix(prev, all_sim_prevs, delta, first_weight)
        ess <- trachomAMIS::calculate_ess(WW)
        cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))
        if(min(ess) >= ESS.R) break
    }

    return(param)
}
