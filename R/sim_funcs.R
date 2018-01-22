#' Dirichlet distributed random numbers
#'
#' Generate n random vectors distributed according to a Dirichlet distribution.
#' Each row of the returned matrix is a random vector that sums to 1.
#'
#' @param n The number of random vectors to generate
#' @param alpha A vector of parameters
#'
#' @return The function returns a matrix with n rows and \code{length(alpha)} columns
#'
#' @source copied from package \code{MCMCpack} to avoid a dependency. That code
#' was taken from Greg's Miscellaneous Functions (gregmisc). His code was based
#' on code posted by Ben Bolker to R-News on 15 Dec 2000.
#'
rdirichlet <- function (n, alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
}

#' Simulate population projection matrices
#'
#' \code{sim_transition} generates a list of simulated population projection matrices from the provided parameters
#' and prior distributions.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stages at start of transition.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param samples The number of matrices to return.
#'
#' @return Always returns a list.
#' @export
#'
sim_transitions <- function(TF, N, P = NULL, priorweight = -1, samples = 1){
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  if(missing(P)){
    # fill in with a uniform prior <- <- <-
    P <- matrix(1/(order+1), nrow=order+1, ncol = order)
  } else {
    if(ncol(P)!=order | nrow(P) != (order+1)) {
      stop("Bad dimensions on P")
    }
  }
  TN <- fill_transitions(TF, N, P, priorweight, returnType = "TN")
  Amats <- list()
  for(i in 1:samples){
    A <- matrix(0, nrow=order, ncol=order)
    for(j in 1:order){
      A[,j] <- rdirichlet(1, TN[,j])[1:order]
    }
    Amats[[i]] <- A + Fmat
  }


  return(Amats)
}