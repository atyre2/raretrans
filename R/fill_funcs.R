#' Combine prior and data for transition matrix
#'
#' \code{fill_transition} returns the expected value of the transition
#' matrix combining observed transitions for one time step and a prior
#'
#' This is just an intermediate version to test package structure. Does
#' not yet use the prior argument P.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param P A matrix of the priors for each column.
#' @param N A vector of observed transitions.
#' @param returnType A character vector describing the desired return value.
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a dirichlet prior
#'   \item T - just the filled in transition matrix
#'   \item TN - the augmented matrix of fates -- use in calculating the CI
#' }
#'
#' @export
#'
fill_transitions <- function(TF, P, N, returnType = "A"){
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  Tfilled <- matrix(NA, nrow=order, ncol=order)
  TN <- matrix(NA, nrow=order+1, ncol = order)
  for (i in 1:order){
    observed <- Tmat[,i] * N[i]
    allfates <- c(observed, N[i]-sum(observed))
    missing <- allfates == 0
    allfates[missing] <- 1
    allfates[!missing] <- allfates[!missing] + sum(missing)
    Tfilled[,i] <- allfates[1:order] / sum(allfates)
    TN[,i] <- allfates
  }
  if (returnType == "A"){
    return(Tfilled + Fmat)
  } else if (returnType == "T") {
    return(Tfilled)
  } else if (returnType == "TN") {
    return(TN)
  } else {
    stop("Bad returntype in fill_transitions()")
  }
}

#' Combine prior and data for fertility matrix
#'
#' \code{fill_fertility} returns the expected value of the fertility
#' matrix combining observed recruits for one time step and a prior
#'
#' This is just an intermediate version to test package structure. Does
#' not yet use the prior argument P.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param P A matrix of the priors for each column.
#' @param N A vector of observed transitions.
#' @param returnType A character vector describing the desired return value.
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a dirichlet prior
#'   \item F - just the filled in fertility matrix
#' }
#'
#' @export
#'
fill_fecundity <- function(TF, P, N, returnType = "A"){
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  Ffilled <- Fmat
  # only change Ffilled if 3,3 is zero and there are adults
  # if no adults, should be 0 for that year.
  if (Fmat[3,3]<0.000001 & N[3] > 0){
    Ffilled[3,3] <- 1 / N[3] # this isn't going to work if there are no adults!
  }

  if (returnType == "A"){
    return(Tmat + Ffilled)
  } else if (returnType == "F"){
    return(Ffilled)
  } else {
    stop("Bad returntype in fill_fecundity()")
  }
}
