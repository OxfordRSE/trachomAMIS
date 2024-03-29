% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AMIS_steps.R
\name{systematic_sample}
\alias{systematic_sample}
\title{Systematic resampling function}
\usage{
systematic_sample(nsamples, weights, log = F)
}
\arguments{
\item{nsamples}{number of samples to draw}

\item{weights}{vector of length equal to the number of particles, containing their weights}

\item{log}{logical indicating if weights are log-weights}
}
\value{
vector of indices of the sampled particles
}
\description{
Implement systematic resampling to reduce variance in weighted particle selection
}
