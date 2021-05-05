#' @export
amis <- function(prevalence_map, transmission_model, n_params, N, IO_file_id, delta = 5, T = 100, target_ess = 250) {
    param<-matrix(NA, ncol=n_params+1+1, nrow=sum(N))  # Matrix for parameter values, + prevalence and weights
    tmp<-rprop0(N[1])    #N[t] random draws of parameters from prior
    x <- tmp[[1]]  # bet
    y <- tmp[[2]]  # constant

    id <- function(t) paste(IO_file_id, sprintf("_it%g", t), sep = "")
    sim_prev <- trachomAMIS::run_transmission_model(seeds = 1:N[1], x, id(1))
    param[1:N[1],1]<-x
    param[1:N[1],2]<-y
    param[1:N[1],3]<-sim_prev

    WW <- trachomAMIS::compute_weight_matrix(prev, sim_prev, delta, first_weight = rep(1, N[1]))
    ess <- trachomAMIS::calculate_ess(WW)
    cat( min(ess),  "", max(ess), "\n")

    GG<-c(0)
    Sigma <- list(NA, 10*T)
    Mean<-list(NA, 10*T)
    PP<-list(NA,T)
    components <- list(GG = GG,
                       Sigma = Sigma,
                       Mean = Mean,
                       PP = PP)
                                        # Set distribution for proposal: Student's t distribution
    proposal=mvtComp(df=3); mixture=mclustMix();
    dprop <- proposal$d
    rprop <- proposal$r

    set.seed(iscen)
    rows_range <- function(t) (sum(N[1:(t-1)])+1):sum(N[1:(t)])
    for (t in 2:T) {
        WW <- update_according_to_ess_value(WW, ess, target_ess)

        parameters <- param[1:sum(N[1:(t-1)]),1:2]
        clustMix <- trachomAMIS::evaluate_mixture(parameters, N[t], WW, mixture)
        param[rows_range(t),1:2] <- trachomAMIS::sample_new_parameters(clustMix, N[t], rprop)
        param[rows_range(t),3] <-trachomAMIS::run_transmission_model(seeds = c((t-1)*N[t]+1:t*N[t]), param[rows_range(t),1], id(t))
        components <- trachomAMIS::update_mixture_components(clustMix, components, t)
        first_weight <- trachomAMIS::compute_prior_proposal_ratio(components, t, T, N, beta = param[,1], constant = param[,2], dprop)
        all_sim_prevs<-param[1:sum(N[1:(t)]),3]
        WW <- trachomAMIS::compute_weight_matrix(prev, all_sim_prevs, delta, first_weight)
        ess <- trachomAMIS::calculate_ess(WW)
        cat( c("min(ESS)=", min(ess),  ", max(ESS)=", max(ess), "\n"))
        if(min(ess) >= target_ess) break
    }

    return(param[,1:3])
}
