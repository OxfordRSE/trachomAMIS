# functions for plotting the empirical CDF

Ecdf <- function(Samples, Weights=NULL) # calculate empirical cdf
{
  if (is.null(Weights)) {Weights=rep(1/length(Samples),length(Samples))}
  Matrix = data.frame(Samples, Weights)
  Matrix = Matrix[order(Matrix$Samples),]
  Matrix$Heights=cumsum(Matrix$Weights)
  return(Matrix)
}

plot_two_ecdfs <- function(sample_prev, data_prev, weights){
  #plots both the data and sample ecdf
  par(mfrow=c(1,1))
  weights = weights/sum(weights)
  plot(NA, xlim=c(0,100), ylim=c(0,1),
       xlab="Infection prevalence (%)",
       ylab = "Empirical cumulative distribution")
  sample_ecdf_matrix <- Ecdf(sample_prev, Weights=weights)
  data_ecdf_matrix <- Ecdf(data_prev)
  lines(stepfun(sample_ecdf_matrix$Samples,
                c(0,sample_ecdf_matrix$Heights)),
        col=4) # red
  lines(stepfun(data_ecdf_matrix$Samples,
                c(0,data_ecdf_matrix$Heights)),
        col=2) # blue
  legend(x=70, y=.1, legend = c("Sampled CDF", "CDF from data"), fill = c(4,2))
}

make_plots_for_ius <- function(paramw, data_prev, plot_folder, IUlist, scenid, grpid){
  for(i in 1:nrow(data_prev)){
    pdf(file=paste0(plot_folder, "/ecdf.", scenid, IUlist[i], ".pdf"))
    par(mfrow=c(1,1))
    sample_prev <- paramw[,"sim_prev"]
    weights <- paramw[,3+i]
    plot_two_ecdfs(sample_prev, data_prev[i,], weights)
    dev.off()
  }
}
