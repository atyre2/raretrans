#' Combine prior and data for transition matrix
#'
#' \code{fill_transition} returns the expected value of the transition
#' matrix combining observed transitions for one time step and a prior
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed transitions.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param returnType A character vector describing the desired return value. Defaults to "T" the transition matrix.
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a dirichlet prior
#'   \item T - just the filled in transition matrix
#'   \item TN - the augmented matrix of fates -- use in calculating the CI or for simulation
#' }
#'
#' @export
#'
fill_transitions <- function(TF, N, P = NULL, priorweight = -1, returnType = "T"){
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
  Tfilled <- matrix(NA, nrow=order, ncol=order)
  TN <- matrix(NA, nrow=order+1, ncol = order)
  for (i in 1:order){
    observed <- Tmat[,i] * N[i]
    if (priorweight > 0 & N[i] > 0){
      P[,i] <- P[,i]*priorweight*N[i]
    }
    allfates <- c(observed, N[i]-sum(observed)) + P[,i]
    # missing <- allfates == 0
    # allfates[missing] <- 1
    # allfates[!missing] <- allfates[!missing] + sum(missing)
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
#' matrix combining observed recruits for one time step and a Gamma prior for each column.
#'
#' Assumes that only one stage reproduces ... needs generalizing.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stage distribution.
#' @param alpha A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param beta A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param returnType A character vector describing the desired return value. Defaults to "F" the fertility matrix
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a Gamma prior
#'   \item F - just the filled in fertility matrix
#'   \item ab - the posterior parameters alpha and beta as a list.
#' }
#'
#' @export
#'
fill_fertility <- function(TF, N, alpha = 0.00001, beta = 0.00001, priorweight = -1, returnType = "F"){
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  if (length(N) != order | sum(is.na(N)) > 0){
    stop("N isn't the correct length or has missing values.")
  }
  if (!(is.numeric(alpha)&is.numeric(beta))){
    stop("alpha or beta must be numeric matrices or single values.")
  }
  # if (sum(is.na(alpha))>0 | sum(is.na(beta)) > 0) {
  #   stop("alpha or beta cannot have missing values")
  # }
  if (length(alpha) != order^2 | length(beta) != order^2) {
    warning("length(alpha | beta) != order^2: only using first value of alpha and beta")
    alpha <- matrix(rep(alpha[1], order^2), nrow = order, ncol = order)
    beta <- matrix(rep(beta[1], order^2), nrow = order, ncol = order)
  }
  Ffilled <- matrix(NA, nrow=order, ncol = order)
  babies_next_year <- sweep(Fmat, 2, N, FUN = "*")
  # test reproducing stages with beta, because of divide by zero issues
  reproducing_stages <- apply(beta, 2, function(x) sum(!is.na(x))) > 0
  # matrix multiplication doesn't preserve the column structure

  if ((all(N[reproducing_stages] > 0) | sum(beta, na.rm = TRUE) > 0)){
    if (priorweight > 0){
      alpha_post <- sweep(alpha, 2, N, FUN = "*")*priorweight + babies_next_year
      beta_post <- sweep(beta, 2, N, FUN = "*")*priorweight + N
    } else {
      alpha_post <- alpha + babies_next_year
      beta_post <- sweep(beta, 2, N, FUN = "+")
    }
    Ffilled <- alpha_post / beta_post
    Ffilled[is.na(Ffilled)] <- 0
  } else {
    stop("No reproducing stages in N, and no positive values in beta.")
  }

  #shouldn't be any missing values in the output
  if (any(is.na(Ffilled))){
    stop("Missing values in filled fertility matrix: check inputs for missing values.")
  }
  if (returnType == "A"){
    return(Tmat + Ffilled)
  } else if (returnType == "F"){
    return(Ffilled)
  } else if (returnType == "ab"){
    return(list(alpha = alpha_post, beta = beta_post))
  } else {
    stop("Bad returntype in fill_fertility()")
  }
}

#' Helper functions for generating different priors
#'
#' Extract the number of individuals in each stage from a dataframe
#' of transitions.
#'
#' @param transitions a dataframe of observations of individuals in different stages
#' @param stage the name of the variable with the stage information
#' @param sort a vector of stage names in the desired order. Default is the order of levels in stage.
#'
#' @return a vector of the counts of observations in each level of stage.
#' @export
#'

get_state_vector <- function (transitions, stage = NULL,
                              sort = NULL)
{
  if (missing(stage)) {
    stage <- "stage"
  }
  nl <- as.list(1:ncol(transitions))
  names(nl) <- names(transitions)
  stage <- eval(substitute(stage), nl, parent.frame())
  if (is.null(transitions[, stage])) {
    stop("No stage column matching ", stage)
  }
  if (missing(sort)) {
    sort <- levels(transitions[[stage]])
  }
  tf <- table(transitions[[stage]])[sort]
  return(as.vector(tf))
}

