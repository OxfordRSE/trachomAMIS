source('AMIS.R')
source('AMIS_steps.R')
source('clusteringSS.R')
#' Define simple "transmission" model
#'
#' Uses the negative binomial model of worm burden W and overdispersion k.
#' p = 1 - (1+W/k)^-k.
#' See Anderson and May (1992).
#' @param seeds a vector of seeds
#' @param parameters a matrix of sampled parameter vectors
#' @return a vector of prevalences 
transmission_model<-function(seeds,parameters) {
  n_samples<-length(seeds)
  prevalences <- rep(NA,n_samples)
  for (i in 1:n_samples) {
    prevalences[i]<-1 - (1+parameters[i,1]/parameters[i,2])^(-parameters[i,2])
  }
  return(prevalences)
}
L<-3
M<-1000
prevalence_map<-matrix(NA,L,M)
for (l in 1:L) {
  prevalence_map[l,]<-rbeta(M,max(1,l-1),max(1,3-l))
}
rownames(prevalence_map)<-c("Here","There","Everywhere")
rprior <- function(n) {
  params<-matrix(NA,n,2)
  colnames(params)<-c("W","k")
  params[,1]<-rexp(n)
  params[,2]<-rexp(n)
  return(params)
}
dprior <- function(x,log=FALSE) {
  if (log) {
    return(sum(dexp(x,log=T)))
  } else {
    return(prod(dexp(x)))
  }
}
prior<-list(rprior=rprior,dprior=dprior)
amis_params<-default_amis_params()
output<-amis(prevalence_map,transmission_model, prior, amis_params)
