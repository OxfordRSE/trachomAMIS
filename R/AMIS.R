#' @export
read_simulated_prevalence <- function(output_file) {
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
