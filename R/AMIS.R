write_model_input <- function(seeds, parameters, input_file) {
    input_params <- cbind(seeds, parameters)
    colnames(input_params) = c("randomgen", "bet")
    write.csv(input_params, file=input_file, row.names=FALSE)
}

#' @export
run_transmission_model <- function(seeds, parameters, input_file, output_file) {
    write_model_input(seeds, parameters, input_file)
    ## Run model
    res <- read.csv(output_file)
    return(100*res[,dim(res)[2]])
}

#' @export
compute_weight_matrix <- function(prev_data, prev_sim, delta, first_weight) {
    n_IUs <- dim(prev_data)[1]
    weight_mat <- matrix(NA, nrow = n_IUs, ncol = length(prev_sim))

    radon_niko_deriv <- function(idx, prev_data_for_IU) {
        f <- length(which((prev_data_for_IU>prev_sim[idx]-delta/2) & (prev_data_for_IU<=prev_sim[idx]+delta/2)))
        g <- sum(first_weight[which((prev_sim>prev_sim[idx]-delta/2) & (prev_sim<=prev_sim[idx]+delta/2))])

    return(f/g)
    }

    for(i in 1:n_IUs) {
        w<-sapply(1:length(prev_sim), radon_niko_deriv, prev_data[i,])
        w <- w*first_weight
        if(sum(w)>0) w<-w/sum(w)
        weight_mat[i,] <- w
    }
    return(weight_mat)
}

#' @export
calculate_ess <- function(weight_mat) {
    ess_for_IU <- function(weights_for_IU) {
        if(sum(weights_for_IU) == 0) return(0)
        return(
            (sum((weights_for_IU)^2))^(-1)
        )
    }

    return(
        apply(weight_mat, 1, ess_for_IU)
    )
}

dprop0<-function(a,b){
    return(dunif(a, min=0.05, max=0.175)*dunif(b, min=0, max=1))
}

#' @export
update_according_to_ess_value <- function(weight_matrix, ess, target_size) {
    rows_to_nullify <- which(ess >= target_size)
    weight_matrix[rows_to_nullify,] <- 0
    return(weight_matrix)
}

#' @export
evaluate_mixture <- function(parameters, nsamples, weight_matrix, mixture) {
    sampled_idx <- sample(
        1:dim(parameters)[1],
        nsamples,
        prob = colSums(weight_matrix),
        replace = T)
    return(
        mixture(parameters[sampled_idx,])
    )
}

#' @export
sample_new_parameters <- function(clustMix, n_samples) {
    x <- c(); y <- c()
    while(length(x)<n_samples){
        compo <- sample(1:clustMix$G,1,prob=clustMix$alpha)
        x1 <- t(
            rprop(1,clustMix$muHat[compo,], clustMix$SigmaHat[,,compo])
        )
        new.param<-as.numeric(x1)
        if(dprop0(new.param[1],new.param[2])>0){
            x<-c(x, new.param[1])
            y<-c(y, new.param[2])
        }
    }
    return(
        list(beta=x, constant=y)
    )
}

#' @export
compute_prior_proposal_ratio <- function(clustMix, t, T, N, beta, constant) {
    ppt <- clustMix$alpha
    muHatt <- clustMix$muHat
    varHatt <- clustMix$SigmaHat
    G <- clustMix$G
    GG<-list(NA,T)
    GG[[t-1]]<-G
    G1<-0; G2<-G
    if(t>2) {
        G1<-sum(sapply(1:(t-2), function(a) GG[[a]]))
        G2<-sum(sapply(1:(t-1), function(a) GG[[a]]))
    }
    Sigma <- list(NA, 10*T)
    Mean<-list(NA, 10*T)
    PP<-list(NA,T)
    for(i in 1:G){
        Sigma[[i+G1]] <- varHatt[,,i]
        Mean[[i+G1]] <- muHatt[i,]
        PP[[i+G1]]<-ppt[i]   ### scale by number of points
    }

    prop.val <- sapply(1:sum(N[1:t]),function(b)  sum(sapply(1:G2, function(g) PP[[g]] * dprop(c(beta[b], constant[b]),mu= Mean[[g]], Sig=Sigma[[g]]))) + dprop0(beta[b], constant[b]))   ## FIX to be just the proposal density ALSO scale by number of points

    first_weight <- sapply(1:sum(N[1:t]), function(b) dprop0(beta[b], constant[b])/prop.val[b])   # prior/proposal

    return(first_weight)
}
