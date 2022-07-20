#' Produce list containing the default AMIS parameters
#' @params histogram logical indicating whether to use the histogram method.
#' @return list containing the default AMIS parameters
#' @export
default_amis_params <- function(histogram=FALSE) {
  amis_params<-list(delta=0.01,nsamples=500,mixture_samples=1000,df=3,target_ess=500,log=F,max_iters=12)
  if (histogram) {amis_params[["breaks"]]<-0:100/100}
  return(amis_params)
}

#' Compute likelihood for each additional simulation across timepoints
#'
#' Calls evaluate likelihood for each timepoint.
#' @param prevalence_map A list with one entry for each timepoint.
#' Each entry must be a list containing objects \code{data} (an L x M matrix of data);
#' and \code{likelihood} a function taking arguments \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.
#' If a likelihood is not specified then it is assumed that
#' the data consist of samples from a geo-statistical model and empirical methods are used.
#' @param An n x timepoints matrix of prevalences simulated from the transmission model.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations (optional).
#' @return A larger array with the likelihoods of the new simulations joined to the existing array \code{likelihoods}.
compute_likelihood<-function(prevalence_map,simulated_prevalences,amis_params,likelihoods=NULL) {
  n_tims<-length(prevalence_map)
  n_locs <- dim(prevalence_map[[1]]$data)[1]
  n_sims <- dim(simulated_prevalences)[1]
  lik<-array(NA,c(n_tims,n_locs,n_sims)) # this way around to avoid abind -- different to elsewhere
  for (t in 1:n_tims) {
    lik[t,,]<-evaluate_likelihood(prevalence_map[[t]],simulated_prevalences[,t],amis_params)
  }
  if (!is.null(likelihoods)) {lik<-array(c(likelihoods,lik),c(n_tims,n_locs,dim(likelihoods)[3]+n_sims))}
  return(lik)
}

#' Evaluate likelihood for each additional simulation for a single timepoint
#'
#' Implements analytical likelihoods if a likelihood function is available; otherwise histogram or empirical
#' likelihoods are generated based on samples from a geostatistical map.
#' @param prevalence_map A list containing objects \code{data} (an L x M matrix of data);
#' and \code{likelihood} a function taking arguments \code{data} (a matrix of data as above),
#' \code{prevalence} (a matrix of output from the transmission model) and optional logical \code{log}, which returns the vector of (log)-likelihoods.
#' @param prev_sim A vector of simulated prevalences
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @return A locations x simulations matrix of (log-)likelihoods.
evaluate_likelihood<-function(prevalence_map,prev_sim,amis_params) {
  locs<-which(!is.na(prevalence_map$data[,1])) # if there is no data for a location, do not update weights.
  f<-matrix(NA,dim(prevalence_map$data)[1],length(prev_sim))
  if (!is.null(prevalence_map$likelihood)) {
    f[locs,]<-t(prevalence_map$likelihood(prevalence_map$data[locs,,drop=FALSE],prev_sim,amis_params[["log"]])) # likelihood function must be vectorised.
  } else {
    if (is.null(amis_params[["breaks"]])) {
      delta<-amis_params[["delta"]]
      for (i in 1:length(prev_sim)) {
        f[,i]<-rowSums(abs(prevalence_map$data[locs,,drop=FALSE]-prev_sim[i])<=delta/2)/delta
      }
    } else {
      breaks<-amis_params[["breaks"]] # NB top entry in breaks must be strictly larger than the largest possible prevalence.
      L<-length(breaks)
      lwr<-breaks[1:(L-1)]
      upr<-breaks[2:L]
      wdt<-upr-lwr
      for (l in 1:L) {
        wh<-which(prev_sim>=lwr[l] & prev_sim<upr[l])
        if (length(wh)>0) {
          f[locs,wh]<-rowSums(prevalence_map$data[locs,,drop=FALSE]>=lwr[l] & prevalence_map$data[locs,,drop=FALSE]<upr[l])/wdt[l]
        }
      }
    }
    if (amis_params[["log"]]) {f<-log(f)}
  }
  return(f)
}

#' Compute weight matrix across timepoints using appropriate method
#'
#' Wrapper function to select appropriate method to calculate weight matrix.
#' @param likelihoods An array with dimension n_tims,n_locs,n_sims -- ie timepoints x locations x simulations.
#' @param simulated_prevalence An n_sims x n_tims matrix containing the simulated prevalence values for each of the
#'     parameter samples. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}.
#' @param first_weight A vector containing the values for the right hand side of
#'     the weight expression. Should be of the same length as the rows in \code{simulated_prevalence}.
#' @return normalised weight matrix.
compute_weight_matrix <- function(likelihoods, simulated_prevalence, amis_params, first_weight) {
  n_tims <- dim(likelihoods)[1]
  n_locs <- dim(likelihoods)[2]
  n_sims <- dim(likelihoods)[3]
  weight_matrix <- matrix(rep(first_weight,n_locs), nrow = n_sims, ncol = n_locs)
  for (t in 1:n_tims) {
    # Update the weights by the latest likelihood (filtering)
    if (is.null(amis_params[["breaks"]])) {
      weight_matrix <- compute_weight_matrix_empirical(t(likelihoods[t,,]),simulated_prevalence[,t],amis_params,weight_matrix)
    } else {
      weight_matrix <- compute_weight_matrix_histogram(t(likelihoods[t,,]),simulated_prevalence[,t],amis_params,weight_matrix)
    }
  }
  # renormalise weights
  if (amis_params[["log"]]) {
    M<-apply(weight_matrix,2,max)
    wh<-which(M>-Inf)
    weight_matrix[,wh]<-weight_matrix[,wh]-rep(M[wh]+log(colSums(exp(weight_matrix[,wh,drop=FALSE]-rep(M[wh],each=n_sims)))),each=n_sims)
  } else {
    S<-colSums(weight_matrix)
    wh<-which(S>0)
    weight_matrix[,wh]<-weight_matrix[,wh]/rep(S[wh],each=n_sims)
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
#' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
#' NB: transpose of slice of array.
#' @param prev_sim A vector containing the simulated prevalence value for each parameter sample.
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#' @param weight_matrix An n_sims x n_locs matrix containing the current values of the weights.
#' @return An updated weight matrix.
compute_weight_matrix_empirical <- function(likelihoods, prev_sim, amis_params, weight_matrix) {
  delta<-amis_params[["delta"]]
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  new_weights<-matrix(ifelse(amis_params[["log"]],-Inf,0),length(prev_sim),length(likelihoods[1,]))
  for (i in 1:length(prev_sim)) {
    wh<-which(abs(prev_sim-prev_sim[i])<=delta/2)
    g_terms<-weight_matrix[wh,locs,drop=FALSE]
    if (amis_params[["log"]]) {
      M<-apply(g_terms,2,max)
      non_zero_locs<-locs[which(M>-Inf)]
      M<-M[which(M>-Inf)]
      new_weights[i,non_zero_locs]<-weight_matrix[i,non_zero_locs]+likelihoods[i,non_zero_locs]-M-log(colSums(exp(g_terms-rep(M,each=length(wh)))))+log(delta)
    } else {
      g<-colSums(g_terms)
      non_zero_locs<-locs[which(g>0)]
      g<-g[which(g>0)]
      new_weights[i,non_zero_locs]<-weight_matrix[i,non_zero_locs]*likelihoods[i,non_zero_locs]/g*delta
    }
  }
  return(new_weights)
}

#' Compute weight matrix using empirical Radon-Nikodym derivative (with fixed breaks)
#'
#' Compute matrix describing the weights for each parameter sampled, for each
#' location. One row per sample, one column per location.  Each weight
#' is computed based on the empirical Radon-Nikodym derivative, taking into account
#' geostatistical prevalence data for the specific location and the prevalence value
#' computed from the transmission model for the specific parameter sample.
#'
#' @param likelihoods An n_sims x n_locs matrix of (log-)likelihoods
#' NB: transpose of slice of array.
#' @param prev_sim A vector containing the simulated prevalence value for each
#'     parameter sample. (double)
#' @param amis_params A list of parameters, e.g. from \code{\link{default_amis_params}}
#'
#' @param weight_matrix A matrix containing the current values of the weights.
#' @return An updated weight matrix.
compute_weight_matrix_histogram<-function(likelihoods, prev_sim, amis_params, weight_matrix) {
  breaks<-amis_params[["breaks"]] # NB top entry in breaks must be strictly larger than the largest possible prevalence.
  L<-length(breaks)
  lwr<-breaks[1:(L-1)]
  upr<-breaks[2:L]
  wdt<-upr-lwr
  locs<-which(!is.na(likelihoods[1,])) # if there is no data for a location, do not update weights.
  new_weights<-matrix(ifelse(amis_params[["log"]],-Inf,0),length(prev_sim),length(likelihoods[1,]))
  for (l in 1:L) {
    wh<-which(prev_sim>=lwr[l] & prev_sim<upr[l])
    if (length(wh)>0) {
      g_terms<-weight_matrix[wh,locs,drop=FALSE]
      if (amis_params[["log"]]) {
        M<-apply(g_terms,2,max)
        non_zero_locs<-locs[which(M>-Inf)]
        M<-M[which(M>-Inf)]
        new_weights[wh,non_zero_locs]<-weight_matrix[wh,non_zero_locs]+likelihoods[wh,non_zero_locs]-rep(M+log(colSums(exp(g_terms-M))),each=length(wh))+log(wdt[l])
      } else {
        g<-colSums(g_terms)
        non_zero_locs<-locs[which(g>0)]
        g<-g[which(g>0)]
        new_weights[wh,non_zero_locs]<-weight_matrix[wh,non_zero_locs]*likelihoods[wh,non_zero_locs]/rep(g,each=length(wh))*wdt[l]
      }
    }
  }
  return(new_weights)
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
  ess<-rep(0,dim(weight_mat)[2])
  if (log) {
    M<-apply(weight_mat,2,max)
    wh<-which(M>-Inf)
    M<-M[wh]
    ess[wh]<-exp(-2*M)*colSums(exp(2*(weight_mat[,wh,drop=FALSE]-rep(M,each=dim(weight_mat)[1]))))^(-1)
  } else {
    S<-colSums(weight_mat^2)
    wh<-which(S>0)
    ess[wh]<-1/S[wh]
  }
  return(ess)
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
    return(M+log(rowSums(exp(weight_matrix[,active_cols,drop=FALSE]-M))))
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
#' @return A list containing \code{params}, an \code{n_samples} x d matrix containing the sampled parameter values and
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
#' Integrating geostatistical maps and infectious disease transmission models
#' using adaptive multiple importance sampling.
#' Renata Retkute, Panayiota Touloupou, Maria-Gloria Basanez,
#' T. Deirdre Hollingsworth, Simon E.F. Spencer
#' Ann. Appl. Stat. 15 (4) 1980 - 1998, December 2021.
#' DOI: https://doi.org/10.1214/21-AOAS1486
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
  q_terms<-matrix(NA,nrow(param),G)
  for (g in 1:G) {
    if (log) {
      q_terms[,g]<-log(probs[[g]])+mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=T)
    } else {
      q_terms[,g]<-probs[[g]]*mnormt::dmt(param,mean=Mean[[g]],S=Sigma[[g]],df=df,log=F)
    }
  }
  if (log) {
    M<-pmax(apply(q_terms,1,max),prior_density)
    return(prior_density - M - log(rowSums(exp(q_terms-M))+exp(prior_density-M)))
  } else {
    return(prior_density/(rowSums(q_terms)+prior_density))
  }
}
