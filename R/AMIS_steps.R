dprop0<-function(x){
    return(dunif(x[1], min=0.05, max=0.175)*dunif(x[2], min=0, max=1))
}

get_initial_parameters<-function(n){
    init_beta_samples <- runif(n, min=0.05, max=0.175)
    init_constant_samples <- runif(n, min=0, max=1)
    return(
        matrix(c(init_beta_samples, init_constant_samples),
               ncol = 2)
        )
}

write_model_input <- function(seeds, parameters, input_file) {
    input_params <- cbind(seeds, parameters)
    colnames(input_params) = c("randomgen", "bet")
    write.csv(input_params, file=input_file, row.names=FALSE)
}

#' @export
run_transmission_model <- function(seeds, parameters, id) {
    input_file <- paste("files/InputBet_", id, ".csv", sep = "")
    write_model_input(seeds, parameters, input_file)
    ## Run model
    output_file <- paste("output/OutputPrev_", id, ".csv", sep = "")
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
sample_new_parameters <- function(clustMix, n_samples, rprop) {
    x <- c(); y <- c()
    while(length(x)<n_samples){
        compo <- sample(1:clustMix$G,1,prob=clustMix$alpha)
        x1 <- t(
            rprop(1,clustMix$muHat[compo,], clustMix$SigmaHat[,,compo])
        )
        new.param<-as.numeric(x1)
        if(dprop0(new.param)>0){
            x<-c(x, new.param[1])
            y<-c(y, new.param[2])
        }
    }
    return(
        matrix(c(x,y), ncol = 2)
    )
}

update_mixture_components <- function(clustMix, components, t) {
    components$GG[t] <- clustMix$G
    G1<-sum(components$GG[1:(t-1)])
    for(i in 1:clustMix$G){
        components$Sigma[[i+G1]] <- clustMix$SigmaHat[,,i]
        components$Mean[[i+G1]] <- clustMix$muHat[i,]
        components$PP[[i+G1]] <- clustMix$alpha[i]   ### scale by number of points
    }
    return(components)
}


#' @export
compute_prior_proposal_ratio <- function(components, param, dprop) {
    PP <- components$PP
    Sigma <- components$Sigma
    Mean <- components$Mean

    G2<-sum(components$GG)
    prop.val <- sapply(1:nrow(param),function(b)  sum(sapply(1:G2, function(g) PP[[g]] * dprop(param[b,],mu= Mean[[g]], Sig=Sigma[[g]]))) + dprop0(param[b,]))   ## FIX to be just the proposal density ALSO scale by number of points

    first_weight <- sapply(1:nrow(param), function(b) dprop0(param[b,])/prop.val[b])   # prior/proposal

    return(first_weight)
}
