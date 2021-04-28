#' @export
read_simulated_prevalence <- function(output_file) {
    res <- read.csv(output_file)
    return(100*res[,dim(res)[2]])
}
