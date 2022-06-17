# SS note: I think these functions were written by Luca Pozzi poz.luc@gmail.com and taken from the ARAMIS R package.
# https://github.com/mrpozzi/ARAMIS

#' Obtain t-distribution density and RNG
#'
#' Function to output list containing density and random number generator fot t-distribution with specific degrees of freedom.
#' @param df degrees of freedom
#' @return list containing "d" (density function) and "r" (RNG) for t-distribution
mvtComp<-function(df=3){
	list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
		mnormt::dmt(xx,mean=mu,S=Sig,df=df,log=log)
		},"r"=function(n=1,mu=0,Sig=1){
			mnormt::rmt(n,mean=mu,S=Sig,df=df)
			})
}

#' Wrapper function to load mclust package and output a function to fit MVN mixtures
#'
#' Loads R package mclust and produces a function that fits multivariate normal (MVN) mixtures with different numbers of components in the mixture.
#' The mixture with the lowest BIC is outputted. 
#' @param G Vector containing the numbers of components in the mixture to try.#
#' @return A function that fits MVN mixtures. 
mclustMix<-function(G=1:10){
    ## LP?:
    ## Attaching mclust because issue with using
    ## mclust::mvn(data, modelName="VVV",). it should point to
    ## mclust::mvnXXX but this function cannot be found. It might because
    ## mvnXXX is called from eval(expt, parent.frame()) in mvn.
    ## See https://github.com/cran/mclust/blob/598224fc49cf8578ade7190ed73a71a51304267d/R/mclust.R#L4532
	require(mclust)
	if (any(as.numeric(G)<=0)) stop("G must be positive")
	
	function(xx){
		
		clustering <- fitMclust(xx,modelName="VVV",G= G)
		
		G <- clustering$G
		
		if(G==1) clustering$parameters$pro <- 1

		return(list(alpha=clustering$parameters$pro, muHat=t(clustering$parameters$mean), SigmaHat=clustering$parameters$variance$sigma,G=G,cluster=clustering$classification))
		}		
}


#' Function to fit MVN normal mixtures and choose the best number of components using BIC.
#'
#' @param xx A matrix with n rows and p columns, where n is the number of observations and p is the dimension
#' @param modelName string containing the model name to pass to mclust
#' @param G vector containing the possible numbers of components to try
#' @return A list containing the best MVN mixture   
fitMclust<-function(xx,modelName="VVV",G= G){
	
	options(warn=-1)
	# set options for EM algorithm
	control <- mclust::emControl(eps=sqrt(.Machine$double.eps))
	
	n <- nrow(xx) # number of observations
	p <- ncol(xx) # dimension
	
	clustering <-Gout <- BIC <- NA # make empty output variables

	if (G[1] == 1) { 
    if (p==1) { # SS note: I added this section to deal with 1D case seperately
  		clustering <- mclust::mvn(modelName = "X", data = xx)
  		BIC <- bic(modelName="X",loglik=clustering$loglik,n=n,d=p,G=1)
    } else { # Fit MVN (ie mixture with just 1 component) and store results as best yet.
  		clustering <- mclust::mvn(modelName = modelName, data = xx)
  		BIC <- bic(modelName=modelName,loglik=clustering$loglik,n=n,d=p,G=1)
    }
		Gout <- 1
		G <- G[-1]
		}
	
	if (p != 1) {
		if (n > p) {
			hcPairs <- hc(modelName="VVV",data=xx) # fit agglomerative hierarchical MVN mixture using mclust
			}else {
				hcPairs <- hc(modelName="EII",data=xx) # if not enough data, fit simpler hierarchical MVN model using mclust
				}
		}else hcPairs <- NULL
	if (p > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G) # evaluare hierarchical clustering tree at group sizes G
	
	for (g in G) {
		
		if (p > 1 || !is.null(hcPairs)) {
			cl <- clss[, as.character(g)]
			}else { # Else perform fail-safe 1D clustering
        #cl <- .qclass(data[subset], as.numeric(g)) # I'm pretty sure this line won't work; what are data and subset? 
				cl <- .qclass(xx, as.numeric(g)) # SS: My attempt to correct the above line
				}
    # convert cluster mapping into a vector of cluster indices.    
		z <- unmap(cl, groups = 1:max(cl))

    if (n>p) {
		   new <- me(modelName=modelName,data=xx,z=z,control=control) # Run EM algorithm to obtain means and covariance matrices of clusters
    } else {
      new <- me(modelName="EII",data=xx,z=z,control=control)
		}
		if(!is.na(new$loglik)){
			
			BICnew <- bic(modelName=modelName,loglik=new$loglik,n=n,d=p,G=g,equalPro=control$equalPro) # calculate BIC if EM has produced a loglikelihood
			
			if(is.na(BIC)){ # if initial attempt has failed, take new in its place
				clustering <- new
				BIC <- BICnew
				Gout <- g
				}else{
					if(BICnew>BIC){ # If new attempt is better than previous best, replace.
						clustering <- new
						BIC <- BICnew
						Gout <- g
						}
					}
			}
		}
	
	options(warn=0)
	
	return(c(clustering,G=Gout))
		
	}
#' Function to split a 1D data vector in k quantiles of size 1/k.
#' Possibly a copy of an internal method in hclust
#' used as failsafe clustering here
#' @param x the data to split into quantiles
#' @param k the number of groups to split the data into
#' @return a length(x) vector of cluster indeces for 1 to k.
.qclass <- function (x, k) 
{
  q <- quantile(x, seq(from = 0, to = 1, by = 1/k))
  cl <- rep(0, length(x))
  q[1] <- q[1] - 1
  for(i in 1:k) 
     cl[x > q[i] & x <= q[i+1]] <- i
  
  return(cl)
}
