#' Fit an MVN mixture model using mclust
#' 
#' Can handle 1D as well as multivariate clustering. \code{dat} must have \code{nrow} observation and \code{ncol}
#'  dimensions, even if the number of dimensions is 1.
#' Uses BIC to determine the best number of components, up to max.components.
#'
#' @param dat a MATRIX or dataframe containing the observations to cluster.
#' @param max.components A postive integer specifying the maximum number of components to fit in the mixture.
#' @return list containing selected output from \code{mclust}:
#' \describe{
#' \item{\code{G}}{the best number of components G.}
#' \item{\code{probs}}{vector of cluster probabilities (mixing weights).}
#' \item{\code{Mean}}{matrix of cluster means.}
#' \item{\code{Sigma}}{array of cluster covariance matrices.}
#' \item{\code{BIC}}{The BIC of the chosen mixture.}
#' \item{\code{ModelName}}{The model name from the package \code{mclust}.}
#' }
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
    cat(g,em$BIC,"\n")
    if (!is.na(em$BIC) && em$BIC>BIC) {
      clustering<-em
      G<-g
      BIC<-em$BIC
    }
  }
  return(c(clustering,G=G))
}
