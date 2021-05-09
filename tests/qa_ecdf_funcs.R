# functions for plotting the empirical CDF

Ecdf <- function(Samples, Weights=NULL) # calculate empirical cdf
{
  if (is.null(Weights)) {Weights=rep(1/length(Samples),length(Samples))}
  Matrix = data.frame(Samples, Weights)
  Matrix = Matrix[order(Matrix$Samples),]
  Matrix$Heights=cumsum(Matrix$Weights)
  return(Matrix)
}

add.points <- function(M,extra) {
  L<-length(extra)
  O<-data.frame(Samples=extra,Weights=rep(0,L),Heights=rep(NA,L))
  rownames(O)<-paste("e",1:length(extra))
  M<-rbind(M,O)
  M<-M[order(M$Samples),]
  for (i in 1:length(M$Samples)) {
    if (is.na(M$Heights[i])) {
      if (i==1) {
        M$Heights[i]<-0
      } else {
        M$Heights[i]<-M$Heights[i-1]
      }
    }
  }
  return(M)
}

TestsCdfs <- function(s1,w1,s2,w2) {
  # add normalization of weights
  w1 <- w1/sum(w1)
  w2 <- w2/sum(w2)
  Ecdf1<-add.points(Ecdf(s1,w1),s2)
  Ecdf2<-add.points(Ecdf(s2,w2),s1)
  dheight<-Ecdf1$Height-Ecdf2$Height
  dwidth<-diff(Ecdf1$Samples)
  #wh.pos<-which(dheight>0)
  #wh.neg<-which(dheight<0)
  wh.pos<-which(dheight[-length(dheight)]>0)
  wh.neg<-which(dheight[-length(dheight)]<0)
  #if (Ecdf1$Height[length(Ecdf1$Height)]!=1) {stop("Error: cdf1 does not reach 1!")}
  # It is one when you print the value - I guess is rounding error?
  #if (Ecdf2$Height[length(Ecdf2$Height)]!=1) {stop("Error: cdf2 does not reach 1!")}
  A.pos<-sum(dheight[wh.pos]*dwidth[wh.pos])
  A.neg<-sum(dheight[wh.neg]*dwidth[wh.neg])
  Distance<-sum((dheight[wh.pos]^2)*dwidth[wh.pos]) + sum((dheight[wh.neg]^2)*dwidth[wh.neg])
  return(c(A.pos,abs(A.neg),max(abs(dheight), na.rm=T), Distance))
}

qa_plots <- function(paramw, grp_map_data, plot_folder, IUlist, scenid, grpid){
  n.ius <- dim(grp_map_data)[1]
  n.data <- dim(grp_map_data)[2] - 1
  ecdf_sqrd_distance <- c()
  ecdf_max_distance <- c()
  for(i in 1:n.ius){
    data_prev <- as.numeric(grp_map_data[i, 2:(n.data+1)])
    sample_prev <- paramw[,2]
    weights <- paramw[,2+i]
    ecdf_dist <- TestsCdfs(sample_prev, weights, data_prev, rep(1, n.data))
    ecdf_max_distance <- c(ecdf_max_distance, ecdf_dist[3])
    ecdf_sqrd_distance <- c(ecdf_sqrd_distance, ecdf_dist[4])
  }
  return(
      data.frame(IUlist, ecdf_max_distance, ecdf_sqrd_distance)
  )
}
