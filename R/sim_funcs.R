#' Dirichlet distributed random numbers
#'
#' Generate n random vectors distributed according to a Dirichlet distribution.
#' Each row of the returned matrix is a random vector that sums to 1.
#'
#' @param n The number of random vectors to generate
#' @param alpha A vector of parameters
#'
#' @return The function returns a matrix with n rows and \code{length(alpha)} columns
#' @export
#'
#' @source copied from package \code{MCMCpack} to avoid a dependency. That code
#' was taken from Greg's Miscellaneous Functions (gregmisc). His code was based
#' on code posted by Ben Bolker to R-News on 15 Dec 2000.
#'
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(stats::rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}

#' Simulate population projection matrices
#'
#' \code{sim_transition} generates a list of simulated population projection matrices from the provided parameters
#' and prior distributions.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stages at start of transition.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param alpha A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param beta A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param samples The number of matrices to return.
#'
#' @return Always returns a list.
#' @export
#'
sim_transitions <- function(TF, N, P = NULL, alpha = 0.00001, beta = 0.00001, priorweight = -1, samples = 1) {
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  if (missing(P)) {
    # fill in with a uniform prior <- <- <-
    P <- matrix(1 / (order + 1), nrow = order + 1, ncol = order)
  } else {
    if (ncol(P) != order | nrow(P) != (order + 1)) {
      stop("Bad dimensions on P")
    }
  }
  TN <- fill_transitions(TF, N, P, priorweight, returnType = "TN")
  # alpha and beta are checked in fill_fertility
  ab_post <- fill_fertility(TF, N, alpha = alpha, beta = beta, priorweight = priorweight, returnType = "ab")
  alpha <- ab_post$alpha # these will now be square matrices
  beta <- ab_post$beta # these will now be square matrices

  Amats <- list()
  for(s in 1:samples){
    T_ <- matrix(0, nrow=order, ncol=order)
    F_ <- matrix(0, nrow=order, ncol=order)

    for(j in 1:order){ # looping over the column of the matrix
      T_[,j] <- rdirichlet(1, TN[,j])[1:order]
      for (i in 1:order){ # loop over rows as well to generate fertilities
        if (!is.na(alpha[i,j])){
          F_[i,j] <- stats::rgamma(1, shape = alpha[i,j], rate = beta[i,j])
        }
      }

    }
    Amats[[s]] <- T_ + F_
  }


  return(Amats)
}

#' Simulate observations from a given matrix
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N integer. Number of plants to assume in the population.
#'
#' @details Currently the function assumes individuals
#' are allowed to reproduce regardless of whether they survive to the next census or
#' not.
#'
#'
#'
#' @return A data.frame with 3 columns for stage, fate, and fertility. Both stage and fate are returned as ordered factors.
#' @export
#'
#' @examples
#' data(calathea, package = "popbio")
#' TF <- splitA(calathea$pooled, r = 1, c = 3:8)
#' sim_observations(TF, 10)
sim_observations <- function(TF, N = 100){
  raretrans:::check_TF(TF)
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  A <- Tmat + Fmat
  eigen_results <- eigen.analysis(A)
  # augment Tmat with row for dead
  Tmat2 <- matrix(NA, nrow = order+1, ncol = order)
  Tmat2[1:order, ] <- Tmat
  Tmat2[order+1, ] <- 1 - colSums(Tmat)
  # check that all columns add to 1
  if(any(abs(colSums(Tmat2)-1) > 1e-5)) stop("column sums not equal to 1")

  if (!is.null(dimnames(Tmat))){
    stage_names <- c(dimnames(Tmat)[[1]], "dead")
  } else {
    stage_names <- c(paste0("stage:", 1:order), "dead")
  }

  results <- data.frame(stage = rep(NA_character_, N),
                        fate = rep(NA_character_, N),
                        fertility = rep(NA_real_, N))
  Ninit <- rmultinom(1, size = N, prob = eigen_results$stable.stage)
  stages <- list()
  fates <- list()
  fertilities <- list()
  for (o in 1:order){
    if (Ninit[o] > 0){
      stages <- append(stages, rep(stage_names[o], Ninit[o]))
      tfate <- rmultinom(1, size = Ninit[o], prob = Tmat2[,o])
      tfate2 <- rep(stage_names, times = tfate)
      tfate3 <- tfate2[sample(Ninit[o])]
      fates <- append(fates, tfate3)
      fertilities <- append(fertilities, rpois(Ninit[o], lambda = Fmat[1, o]))
    }
  }
  results$stage <- ordered(unlist(stages), stage_names[1:order])
  results$fate <- ordered(unlist(fates), stage_names[1:order])
  results$fertility <- unlist(fertilities)

  return(results)
}
