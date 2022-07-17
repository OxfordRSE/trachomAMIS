dprop0 <- function(x) {
  return(dunif(x[1], min = 0.05, max = 0.175) * dunif(x[2], min = 0, max = 1))
}

#' Write input file required by transmission model
#'
#' Write a 2 columns csv file named INPUT_FILE to disk. The first column SEEDS
#' contains the seed for each parameter value. Second column BETA contains
#' corresponding parameter values.
#'
#' @param seeds A vector containing seed value for each parameter sample
#'     (double)
#' @param beta A vector containing samples of the beta parameter (double)
#' @param input_file The name of the input file (string)
write_model_input <- function(seeds, beta, input_file) {
  input_params <- cbind(seeds, beta)
  colnames(input_params) <- c("randomgen", "bet")
  write.csv(input_params, file = input_file, row.names = FALSE)
}

#' Compute weight matrix
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' Implementation Unit (IU). One row per sample, one column per IU.  Each weight
#' is computed based on the Randon-Nikodim derivative, taking into account
#' geostatistical prevalence data for the specific IU and the prevalence value
#' computed from the transmission model for the specific parameter sample.
#'
#' @param prev_data The geostatistical prevalence data. A n x m matrix where n
#'     is the number of IUs and m the number of prevalence samples. (double)
#' @param prev_sim A vector containing the siumulated prevalence value for each
#'     parameter sample. (double)
#' @param delta A number. (double)
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as \code{prev_sim}.
# SS note: make alternative versions of this function for other ERN derivatives.
compute_weight_matrix <- function(prev_data, prev_sim, amis_params, first_weight) {
  n_IUs <- dim(prev_data)[1]
  weight_mat <- matrix(NA, nrow = n_IUs, ncol = length(prev_sim))
  delta<-amis_params[["delta"]]
  # define function to calculate empirical RN derivative from Touloupou, Retkute, Hollingsworth and Spencer (2020)
  radon_niko_deriv <- function(idx, prev_data_for_IU, log=amis_params[["log"]]) {
    f <- length(which((prev_data_for_IU > prev_sim[idx] - delta / 2) & (prev_data_for_IU <= prev_sim[idx] + delta / 2)))
    g_terms <- first_weight[which((prev_sim > prev_sim[idx] - delta / 2) & (prev_sim <= prev_sim[idx] + delta / 2))]
    if (log) {
      M<-max(g_terms)
      return(log(f)-M-log(sum(exp(g_terms-M))))
    } else {
      return(f/sum(g_terms))
    }
  }
  for (i in 1:n_IUs) {
    w <- sapply(1:length(prev_sim), radon_niko_deriv, prev_data[i, ])
    if (amis_params[["log"]]) {
      w <- w + first_weight
      M<-max(w)
      S<-M+log(sum(exp(w-M)))
      if (S>-Inf) {w <- w - S}
    } else {
      w <- w * first_weight
      if (sum(w) > 0) {w <- w / sum(w)}
      # NB it doesn't matter if all the simulations have weight zero in an early iteration, as long as this is not true for all active IUs.
    }
    weight_mat[i, ] <- w
  }
  return(weight_mat)
}

#' Compute the current effective sample size
#'
#' This function returns the effective sample size (ESS) for each Implementation
#' Unit (IU). For each column of the weight matrix WEIGHT_MAT, the ESS component
#' is computed as
#'
#' \deqn{\left(\sum_{i=1}^{N} w_{i}^2 \right)^{-1}},
#'
#' where \eqn{$N$} is the number of sampled parameter values for each IU.
#'
#' @param weight_mat The weight matrix. A n x m matrix with n the number of IUS
#'     and m the number of sampled parameter values.
#' @param log logical indicating if the weights are on the log scale. 
#' @return A vector containing the ESS value for each IU.
#'
#' @seealso \code{\link{compute_weight_matrix}}
calculate_ess <- function(weight_mat,log) {
  if (log) {
    ess_for_IU <- function(weights_for_IU) {
      M <- max(weights_for_IU)
      if (M == -Inf) {
        return(0)
      } else {
        return((exp(2*M)*sum(exp(2*(weights_for_IU-M))))^(-1))
      }
    }
  } else {
    ess_for_IU <- function(weights_for_IU) {
      if (sum(weights_for_IU) == 0) {
        return(0)
      } else {
        return(sum(weights_for_IU^2)^(-1))
      }
    }
  }
  return(apply(weight_mat, 1, ess_for_IU))
}
#' Calculate sum of weight matrix for active locations
#'
#' This function sums the rows of the weight matrix WEIGHT_MATRIX for which
#' the effective sample size ESS is below a target size TARGET_SIZE.
#'
#' @param weight_matrix The weight_matrix as returned by
#'     \link{compute_weight_matrix}
#' @param ess The effective sample size vector as returned by
#'     \link{calculate_ess}
#' @param target_size A number representing the target size for the sample.
#' @param log A logical indicating if the weights are logged.
#' @return The updated weight matrix. Its size is unchanged.
update_according_to_ess_value <- function(weight_matrix, ess, target_size,log) {
  active_rows <- which(ess < target_size)
  if (log) {
    M<-apply(weight_matrix[active_rows,,drop=FALSE],2,max)
    return(M+log(colSums(exp(weight_matrix[active_rows,,drop=FALSE]-rep(M,each=length(active_rows))))))
  } else {
    return(colSums(weight_matrix[active_rows,,drop=FALSE]))
  }
}
#' Systematic resampling function
#' Implement systematic resampling to reduce variance in weighted particle selection
#' @param nsamples number of samples to draw
#' @param weights vector of length equal to the number of particles, containing their weights
#' @param log logical indicating if weights are log-weights
#' @return vector of indices of the sampled particles 
systematic_sample <- function(nsamples,weights,log=F) {
  if (log) {
    M<-max(weights)
    log_sum_weights<-M+log(sum(exp(weights-M)))
    cum <- cumsum(exp(weights-log_sum_weights))
  } else {
    cum <- cumsum(weights)/sum(weights) # cumulative sum of normalised weights
  }
  u <- runif(1)/nsamples+0:(nsamples-1)/nsamples
  return(1+matrix(rep(u,length(weights))>rep(cum,each=nsamples),nsamples,length(weights))%*%matrix(1,length(weights),1))
}
#' Fit mixture to weighted sample
#' Weights are implemented by using systematic resampling to obtain an unweighted set of parameters
#' An unweighted mixture is then fitted using \code{fit_mixture}.
#' @param parameters A matrix m x d matrix containing the sampled values for the
#'     d parameters.
#' @param nsamples The number of parameter to resample as data to fit the mixture to.
#' @param weights A vector of weights
#' @param log logical indicating if weights are logged.
#' @return A list of the mixture components (see function \code{\link{fit_mixture}})
#'     \describe{
#'       \item{\code{alpha}}{The mixture weights}
#'       \item{\code{muHat}}{The means of the components}
#'       \item{\code{SigmaHat}}{The covariance matrices of the components}
#'       \item{\code{G}}{Number of components}
#'       \item{\code{cluster}}{clustering IDs}
#'     }
#'
#' @seealso \code{\link{fit_mixture}}
weighted_mixture <- function(parameters, nsamples, weights, log=F) {
  sampled_idx <- systematic_sample(nsamples,weights,log)
  return(fit_mixture(parameters[sampled_idx, ]))
}

#' Sample new parameters
#'
#' This function generates NSAMPLES new model parameter values according the
#' t distribution with \code{df} degrees of freedom and the mixture components mixture.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{evaluate_mixture}}
#' @param n_samples A number of new parameters to sample (integer)
#' @param df The degrees of freedom for the t-dsitributed proposal distribution.
#' @param prior list containing the functions rprior and dprior
#' @param log A logical indicating if densities 
#' @return A list containing params, A NSAMPLES x d matrix containing the sampled parameter values and
#' prior_density, the corresponding vector of prior densities.
#'
#' @seealso \code{\link{fit_mixture}}
sample_new_parameters <- function(mixture, n_samples, df=3, prior, log) {
  prior_density<-rep(NA,n_samples) 
  i<-1
  while (i <= n_samples) {
    compo <- sample(1:mixture$G, 1, prob = mixture$probs)
    proposal <- mnormt::rmt(1,mean=mixture$Mean[,compo],S=mixture$Sigma[,,compo],df=df)
    density_of_proposal <- prior$dprior(proposal,log=log)
    if (!is.na(proposal) && ((log && density_of_proposal>-Inf) || (!log && density_of_proposal>0))) {
      if (i==1) {params<-matrix(NA,n_samples,length(proposal))}
      params[i,]<-proposal
      prior_density[i]<-density_of_proposal
      i<-i+1
    }
  }
  return(list(params=params,prior_density=prior_density))
}

#' Update the components of the mixture
#'
#' This function updates the components of the mixture COMPONENTS according to
#' the current mixture CLUSTMIX generated at iteration T.
#'
#' @param clustMix A list of mixture components as returned by
#'     \code{\link{evaluate_mixure}}
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{GG}}{A numeric vector}
#'       \item{\code{Sigma}}{A list}
#'       \item{\code{Mean}}{A list}
#'       \item{\code{PP}}{A list}
#'     }
#' @param t The current iteration index (integer)
#' @return The updated \code{components} list
#'
#' @seealso \code{\link{evaluate_mixture}}, \code{\link{mclustMix}}
update_mixture_components <- function(clustMix, components, t) {
  components$GG[t] <- clustMix$G
  G1 <- sum(components$GG[1:(t - 1)])
  for (i in 1:clustMix$G) {
    components$Sigma[[i + G1]] <- clustMix$SigmaHat[, , i]
    components$Mean[[i + G1]] <- clustMix$muHat[i, ]
    components$PP[[i + G1]] <- clustMix$alpha[i] ### scale by number of points
  }
  return(components)
}


#' Compute the prior/proposal ratio
#'
#' This function returns the ratio between the prior and proposal distribution
#' for each sampled parameter value (i.e. each row in PARAM).
#'
#' This function returns
#' \deqn(\frac{\pi(\theta _j)}{\frac{1}{N_1 + ... + N_t}\sum_{u=1}^{t}N_u \phi_u(\theta _j)})
#' See step (4) if the AMIS algorithm in
#'
#' Integrating Geostatistical Maps And Transmission Models Using Adaptive
#' Multiple Importance Sampling
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' medRxiv 2020.08.03.20146241;
#' doi: \url{https://doi.org/10.1101/2020.08.03.20146241}
#'
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{GG}}{A numeric vector}
#'       \item{\code{Sigma}}{A list}
#'       \item{\code{Mean}}{A list}
#'       \item{\code{PP}}{A list}
#'     }
#' @param param A 2 column matrix containing the values for the two parameters.
#' @param dprop The proposal distribution as returned by \code{mvtComp()$d}.
#' @return A vector contaning the prior/proposal ratio for each row in
#'     \code{param}
compute_prior_proposal_ratio <- function(components, param, dprop) {
  PP <- components$PP
  Sigma <- components$Sigma
  Mean <- components$Mean

  G2 <- sum(components$GG)
  prop.val <- sapply(1:nrow(param), function(b) sum(sapply(1:G2, function(g) PP[[g]] * dprop(param[b, ], mu = Mean[[g]], Sig = Sigma[[g]]))) + dprop0(param[b, ])) ## FIX to be just the proposal density ALSO scale by number of points

  first_weight <- sapply(1:nrow(param), function(b) dprop0(param[b, ]) / prop.val[b]) # prior/proposal

  return(first_weight)
}
