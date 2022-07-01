#' Compute weight matrix using appropriate method
#' 
#' Wrapper function to select appropriate method to calculate weight matrix.
#' @param prevalence_map A list with one entry for each timepoint.
#' Each entry is a list containing objects \code{data} (an L x M matrix of data); \code{likelihood} a function taking a row of data and the output from the transmission
#' model as arguments (and logical \code{log}) and returning the (log)-likelihood. If a likelihood is not specified then it is assumed that
#' the data consist of samples from a geo-statistical model and empirical methods are used.
#' @param simulated_prevalence An N x timepoints matrix containing the simulated prevalence values for each of the
#'     N parameter samples. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as the rows in \code{simulated_prevalence}.
compute_weight_matrix <- function(prevalence_map, simulated_prevalence, amis_params, first_weight) {
  timepoints<-length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  n_sims <- dim(simulated_prevalence)[1]
  weight_matrix <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs)
  for (t in 1:timepoints) {
    # Update the weights by the latest likelihood (filtering)
    if (is.null(amis_params[["breaks"]])) {
      weight_matrix <- compute_weight_empirical(prevalence_map[[t]],simulated_prevalence[,t],amis_params,weight_matrix)
    } else {
      weight_matrix <- compute_weight_histogram(prevalence_map[[t]],simulated_prevalence[,t],amis_params,weight_matrix)
    }
  }
  # renormalise weights
  if (amis_params[["log"]]) {
    M<-apply(weight_matrix,2,max)
    wh<-which(M>-Inf)
    weight_matrix[,wh]<-weight_matrix[,wh]-M[wh]-log(colSums(exp(weight_matrix[,wh]-rep(M,each=n_sims))))
  } else {
    S<-colSums(weight_matrix)
    wh<-which(S>0)
    weight_matrix[,wh]<-weight_matrix[,wh]/S[wh]
  }
  return(weight_matrix)
}

#' Compute weight matrix using empirical Radon-Nikodym derivative
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based on the empirical Radon-Nikodym derivative, taking into account
#' geostatistical prevalence data for the specific location and the prevalence values
#' computed from the transmission model for the specific parameter sample.
#'
#' @param prevalence_map The geostatistical prevalence data for a single timepoint. 
#' A list containing \code{data}, an L x M matrix where L is the number of locations
#'  and M the number of prevalence samples; and \code{likelihood}, an optional likelihood function.
#' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#' @param weight_matrix A matrix containing the current values of the weights.
#' @return An updated weight matrix.
compute_weight_matrix_empirical <- function(prevalence_map, prev_sim, amis_params, weight_matrix) {
  delta<-amis_params[["delta"]]
  # define function to calculate empirical RN derivative from Touloupou, Retkute, Hollingsworth and Spencer (2020)
  radon_niko_deriv <- function(idx, prev_data_for_IU, weight_vector, log=amis_params[["log"]]) {
    if (is.null(prevalence_map$likelihood)) {
      f <- length(which((prev_data_for_IU >= prev_sim[idx] - delta) & (prev_data_for_IU <= prev_sim[idx] + delta)))
    } else {
      f <- likelihood(prev_data_for_IU, prev_sim[idx],log)
    }
    g_terms <- weight_vector[which((prev_sim >= prev_sim[idx] - delta) & (prev_sim <= prev_sim[idx] + delta))]
    if (log) {
      M<-max(g_terms)
      return(weight_vector[idx]+log(f)-M-log(sum(exp(g_terms-M))))
    } else {
      return(weight_vector[idx]*f/sum(g_terms))
    }
  }
  locs<-which(!is.na(prevalence_map$data[,1])) # if there is no data for a location, do not update weights.
  for (i in locs) {
    weight_matrix[,i] <- sapply(1:length(prev_sim), radon_niko_deriv, prev_data_for_IU=prevalence_map$data[i, ], weight_vector=weight_matrix[,i])
  }
  return(weight_matrix)
}

#' Compute weight matrix using empirical Radon-Nikodym derivative (with fixed breaks)
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based on the empirical Radon-Nikodym derivative, taking into account
#' geostatistical prevalence data for the specific location and the prevalence value
#' computed from the transmission model for the specific parameter sample.
#'
#' @param prevalence_map The geostatistical prevalence data. An L x M matrix where L
#'     is the number of locations and M the number of prevalence samples. (double)
#' @param prev_sim A vector containing the simulated prevalence value for each
#'     parameter sample. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#'
#' @param weight_matrix A matrix containing the current values of the weights.
#' @return An updated weight matrix.
compute_weight_matrix_histogram<-function(prevalence_map, prev_sim, amis_params, weight_matrix) {
  breaks<-amis_params[["breaks"]] # NB top entry in breaks must be strictly larger than the largest possible prevalence. 
  if (min(breaks)>0) {breaks<-c(0,breaks)}
  L<-length(breaks)
  lwr<-breaks[1:(L-1)]
  upr<-breaks[2:L]
  wdt<-upr-lwr
  locs<-which(!is.na(prevalence_map$data[,1])) # if there is no data for a location, do not update weights.
  for (l in 1:L) {
    wh<-which(prev_sim>=lwr[l] & prev_sim<upr[l])
    if (length(wh)>0) {
      for (i in locs) {
        if (is.null(prevalence_map$likelihood)) {
          f<-length(which(prevalence_map$data[i,]>=lwr[l] & prevalence_map$data[i,]<upr[l]))/wdt[l]
        } else {
          like<-function(idx,prev_data,prev_sim) {return(likelihood(prev_data[idx,],prev_sim[idx],log=amis_params[["log"]]))}
          f<-sapply(wh,like,prev_data=prevalence_map$data,prev_sim=prev_sim)
        }
        g_terms<-weight_matrix[wh,i]
        if (amis_params[["log"]]) {
          M<-max(g_terms)
          weight_matrix[wh,i]<-weight_matrix[wh,i]+log(f)-M-log(sum(exp(g_terms-M)))+log(wdt[l])
        } else {
          weight_matrix[wh,i]<-weight_matrix[wh,i]*f/sum(g_terms)*wdt[l]
        }
      }
    }
  }
  return(weight_matrix)
}

#' Compute the current effective sample size
#'
#' This function returns the effective sample size (ESS) for each location. 
#'  For each column of the weight matrix \code{weight_mat}, the ESS component is computed as
#' \deqn{(\sum_{i=1}^{N} w_{i}^2 )^{-1}}
#' where \eqn{N} is the number of sampled parameter values for each location.
#'
#' @param weight_mat The weight matrix. A N x L matrix with L the number of locations
#'     and N the number of sampled parameter values.
#' @param log logical indicating if the weights are on the log scale. 
#' @return A vector containing the ESS value for each location.
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
  return(apply(weight_mat, 2, ess_for_IU))
}
#' Calculate sum of weight matrix for active locations
#'
#' This function sums the rows of the weight matrix \code{weight_matrix} for which
#' the effective sample size ESS is below a target size \code{target_size}.
#'
#' @param weight_matrix The weight_matrix as returned by
#'     \link{compute_weight_matrix}
#' @param ess The effective sample size vector as returned by
#'     \link{calculate_ess}
#' @param target_size A number representing the target size for the sample.
#' @param log A logical indicating if the weights are logged.
#' @return Vector containing the row sums of the active columns of the weight matrix.
update_according_to_ess_value <- function(weight_matrix, ess, target_size,log) {
  active_cols <- which(ess < target_size)
  if (log) {
    M<-apply(weight_matrix[,active_cols,drop=FALSE],1,max)
    wh<-which(M==-Inf)
    M[wh]<-0
    return(M+log(rowSums(exp(weight_matrix[,active_cols,drop=FALSE]-rep(M,length(active_cols))))))
  } else {
    return(rowSums(weight_matrix[,active_cols,drop=FALSE]))
  }
}
#' Systematic resampling function
#' 
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
#' 
#' Weights are implemented by using systematic resampling to obtain an unweighted set of parameters.
#' An unweighted mixture is then fitted using \code{fit_mixture}.
#' @param parameters An N x d matrix containing the sampled values for the d parameters.
#' @param nsamples The number of parameter to resample as data to fit the mixture to.
#' @param weights A vector of weights with length N.
#' @param log logical indicating if weights are logged.
#' @return A list of the mixture components (see function \code{\link{fit_mixture}})
#'     \describe{
#'       \item{\code{probs}}{The mixture weights}
#'       \item{\code{Mean}}{The means of the components}
#'       \item{\code{Sigma}}{The covariance matrices of the components}
#'       \item{\code{G}}{Number of components}
#'       \item{\code{BIC}}{BIC of fitted mixture}
#'       \item{\code{ModelName}}{Model name from package mclust}
#'     }
#'
#' @seealso \code{\link{fit_mixture}}
weighted_mixture <- function(parameters, nsamples, weights, log=F) {
  sampled_idx <- systematic_sample(nsamples,weights,log)
  return(fit_mixture(parameters[sampled_idx,,drop=FALSE]))
}

#' Sample new parameters
#'
#' This function generates \code{nsamples} new model parameter values according the
#' t-distribution with \code{df} degrees of freedom and the mixture components \code{mixture}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{evaluate_mixture}}
#' @param n_samples A number of new parameters to sample (integer)
#' @param df The degrees of freedom for the t-dsitributed proposal distribution.
#' @param prior list containing the functions \code{rprior} and \code{dprior}
#' @param log A logical indicating if densities 
#' @return A list containing \code{params}, an \code{nsamples} x d matrix containing the sampled parameter values and
#' \code{prior_density}, the corresponding vector of prior densities.
#'
#' @seealso \code{\link{fit_mixture}}
sample_new_parameters <- function(mixture, n_samples, df, prior, log) {
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
#' This function updates the mixture \code{components} according to
#' the current mixture \code{mixture} generated at iteration \code{t}.
#'
#' @param mixture A list of mixture components as returned by
#'     \code{\link{fit_mixture}}
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{G}}{A numeric vector containing the number of components from each AMIS iteration}
#'       \item{\code{Sigma}}{A list of covariance matrices for each component}
#'       \item{\code{Mean}}{A list of means for each component}
#'       \item{\code{probs}}{A list probability weights for each component}
#'     }
#' @param t The current iteration index (integer)
#' @return The updated \code{components} list
#'
#' @seealso \code{\link{evaluate_mixture}}, \code{\link{fit_mixture}}
update_mixture_components <- function(mixture, components, t) {
  components$G[t] <- mixture$G
  G_previous <- sum(components$G[1:(t - 1)]) # Number of pre-existing components
  for (i in 1:mixture$G) {
    components$Sigma[[i + G_previous]] <- mixture$Sigma[, , i]
    components$Mean[[i + G_previous]] <- mixture$Mean[,i]
    components$probs[[i + G_previous]] <- mixture$probs[i] ### scale by number of points if nsamples varies by iteration
  }
  return(components)
}
#' Compute the prior/proposal ratio
#'
#' This function returns the ratio between the prior and proposal distribution
#' for each sampled parameter value (i.e. each row in \code{param}).
#' This function returns the first weight
#' See step (4) of the AMIS algorithm in
#' Integrating Geostatistical Maps And Transmission Models Using Adaptive
#' Multiple Importance Sampling
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' medRxiv 2020.08.03.20146241;
#' doi: https://doi.org/10.1101/2020.08.03.20146241
#'
#' @param components A list of mixture components made of
#'     \describe{
#'       \item{\code{G}}{A numeric vector containing the number of components from each AMIS iteration}
#'       \item{\code{Sigma}}{A list of covariace matrices from each component}
#'       \item{\code{Mean}}{A list of means from each component}
#'       \item{\code{probs}}{A list of probability weights for each component (unnormalised)}
#'     }
#' @param param A matrix containing the sampled parameter vectors.
#' @param prior_density Vector containing the prior density of each sampled parameter vector.
#' @param df The degrees of freedom for the t-distributed proposal distribution.
#' @param log A logical indicating whether to work on the log scale.
#' @return A vector containing the prior/proposal ratio for each row in
#'     \code{param}
compute_prior_proposal_ratio <- function(components, param, prior_density, df, log) {
  probs <- components$probs # /sum(unlist(components$probs)) # to normalise?
  Sigma <- components$Sigma
  Mean <- components$Mean
  G <- sum(components$G)
  if (log) {
    prop.val <- sapply(1:nrow(param), function(b) {
      component_densities <- sapply(1:G, function(g) {log(probs[[g]])+mnormt::dmt(param[b, ], mean = Mean[[g]], S = Sigma[[g]],df=df,log=T)})
      M <- max(component_densities,prior_density[b])
      return(M+log(sum(exp(component_densities-M))+exp(prior_density[b]-M)))
    }) # Assumes number of points is equal in each iteration.
    first_weight <- prior_density - prop.val # prior/proposal
  } else {
    prop.val <- sapply(1:nrow(param), function(b) {sum(sapply(1:G, function(g) {probs[[g]]*mnormt::dmt(param[b, ], mean = Mean[[g]], S = Sigma[[g]],df=df)})) + prior_density[b]}) # Assumes number of points is equal in each iteration.
    first_weight <- prior_density/prop.val # prior/proposal
  }
  return(first_weight)
}
