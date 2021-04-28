
mvtComp<-function(df=3){
	list("d"=function(xx,mu=rep(0,ncol(xx)),Sig=diag(1,ncol(xx),ncol(xx)),log=FALSE){
		dmt(xx,mean=mu,S=Sig,df=df,log=log)
		},"r"=function(n=1,mu=0,Sig=1){
			rmt(n,mean=mu,S=Sig,df=df)
			})
}

mclustMix<-function(G=1:10){
	
	if (any(as.numeric(G)<=0)) stop("G must be positive")
	if(!require("mclust")) stop("Package mclust is not installed...")
	
	function(xx){
		
		clustering <- fitMclust(xx,modelName="VVV",G= G)
		
		G <- clustering$G
		
		if(G==1) clustering$parameters$pro <- 1

		return(list(alpha=clustering$parameters$pro, muHat=t(clustering$parameters$mean), SigmaHat=clustering$parameters$variance$sigma,G=G,cluster=clustering$classification))
		}		
}



fitMclust<-function(xx,modelName="VVV",G= G){
	
	options(warn=-1)
	
	control <- emControl(eps=sqrt(.Machine$double.eps))
	
	n <- nrow(xx)
	p <- ncol(xx)
	
	clustering <-Gout <- BIC <- NA
	
	if (G[1] == 1) {
		clustering <- mvn(modelName = modelName, data = xx)
		BIC <- bic(modelName=modelName,loglik=clustering$loglik,n=n,d=p,G=1)
		Gout <- 1
		G <- G[-1]
		}
	
	if (p != 1) {
		if (n > p) {
			hcPairs <- hc(modelName="VVV",data=xx)
			}else {
				hcPairs <- hc(modelName="EII",data=xx)
				}
		}else hcPairs <- NULL
	if (p > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
	
	for (g in G) {
		
		if (p > 1 || !is.null(hcPairs)) {
			cl <- clss[, as.character(g)]
			}else {
				cl <- .qclass(data[subset], as.numeric(g))
				}
		z <- unmap(cl, groups = 1:max(cl))

		new <- me(modelName=modelName,data=xx,z=z,control=control)
		
		if(!is.na(new$loglik)){
			
			BICnew <- bic(modelName=modelName,loglik=new$loglik,n=n,d=p,G=g,equalPro=control$equalPro)
			
			if(is.na(BIC)){
				clustering <- new
				BIC <- BICnew
				Gout <- g
				}else{
					if(BICnew>BIC){
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


