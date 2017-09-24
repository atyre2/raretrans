#' Combine prior and data for transition matrix
#'
#' \code{fill_transition} returns the expected value of the transition
#' matrix combining observed transitions for one time step and a prior
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed transitions.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param returnType A character vector describing the desired return value.
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
fill_transitions <- function(TF, N, P = NULL, priorweight = -1, returnType = "A"){
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

## helper functions for generating different priors
#' Extract the number of individuals in each stage from a dataframe
#' of transitions
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

## testing
# df_5 %>%
#   mutate(stage = factor(stage, levels=c("p","j","a","m")),
#          fate = factor(next_stage, levels=c("p","j","a","m"))) %>%
#   filter(POPNUM==209, year==1) %>%
#   as.data.frame() %>%
#   get_state_vector(sort=c("p","j","a"))
