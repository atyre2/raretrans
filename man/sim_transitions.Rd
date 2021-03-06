% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_funcs.R
\name{sim_transitions}
\alias{sim_transitions}
\title{Simulate population projection matrices}
\usage{
sim_transitions(
  TF,
  N,
  P = NULL,
  alpha = 1e-05,
  beta = 1e-05,
  priorweight = -1,
  samples = 1
)
}
\arguments{
\item{TF}{A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.}

\item{N}{A vector of observed stages at start of transition.}

\item{P}{A matrix of the priors for each column. Defaults to uniform.}

\item{alpha}{A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.}

\item{beta}{A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.}

\item{priorweight}{total weight for each column of prior as a percentage of sample size or 1 if negative}

\item{samples}{The number of matrices to return.}
}
\value{
Always returns a list.
}
\description{
\code{sim_transition} generates a list of simulated population projection matrices from the provided parameters
and prior distributions.
}
