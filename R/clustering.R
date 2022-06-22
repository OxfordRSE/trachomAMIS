#' Fit an MVN mixture model using mclust
#' Can handle 1D as well as multivariate clustering. dat must have nrow observation and ncol dimensions, even if the number of dimensions is 1.
#' Uses BIC to determine the best number of components, up to max.components.
#'
#' @param dat a MATRIX or dataframe containing the observations to cluster.
#' @param max.components A postive integer specifying the maximum number of components to fit in the mixture.
#' @return list containing the mclust output and the best number of components G.
fit_mixture<-function(dat,max.components=10) {
  require(mclust)
  n<-nrow(dat)
  d<-ncol(dat)
  if (n<d) {stop("Not enough observations to fit mixture model.\n")}
  max.components<-min(n,max.components-1) # or even smaller?
  # Start by fitting one group
  G<-1 # number of groups
  if (d==1) {
    modelName<-"X" 
  } else {
    modelName<-"XXX"
  }
  clustering<-mclust::mvn(modelName=modelName,data=dat)
  BIC <- bic(modelName="X",loglik=clustering$loglik,n=n,d=d,G=1)
  # fit agglomerative clustering model
  if (d==1) {
    modelName<-"V"
  } else {
    modelName <- "VVV"
  }
  hcPairs <- hc(modelName=modelName,data=dat)
  cut.tree <- hclass(hcPairs,2:max.components)
  for (g in 2:max.components) {
    z<-unmap(cut.tree[,g-1]) # extract cluster indices
    # Run EM algorithm
    em <- me(modelName,dat,z)
    em$BIC <- bic(modelName,em$loglik,n,d,g)
    if (!is.na(em$BIC) && em$BIC>BIC) {
      clustering<-em
      G<-g
      BIC<-em$BIC
    }
  }
  return(list(G=G,probs=clustering$parameters$pro,Mean=clustering$parameters$mean,Sigma=clustering$parameters$variance$sigma,
    BIC=BIC,modelName=clustering$modelName))
}
