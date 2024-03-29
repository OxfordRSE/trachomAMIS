% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AMIS_steps.R
\name{update_according_to_ess_value}
\alias{update_according_to_ess_value}
\title{Calculate sum of weight matrix for active locations}
\usage{
update_according_to_ess_value(weight_matrix, ess, target_size, log)
}
\arguments{
\item{weight_matrix}{The weight_matrix as returned by
\link{compute_weight_matrix}}

\item{ess}{The effective sample size vector as returned by
\link{calculate_ess}}

\item{target_size}{A number representing the target size for the sample.}

\item{log}{A logical indicating if the weights are logged.}
}
\value{
Vector containing the row sums of the active columns of the weight matrix.
}
\description{
This function sums the rows of the weight matrix \code{weight_matrix} for which
the effective sample size ESS is below a target size \code{target_size}.
}
