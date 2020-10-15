#' Construct projection matrix models using transition frequency tables
#'
#' Construct an age or stage-structure projection model from a transition table listing stage in time t, fate in time t+1, and one or more individual fertility columns.
#' This version is modified to allow for missing transitions in the data.
#'
#' @param transitions a stage-fate data frame with stage or age class in the current census, fate in the subsequent census, and one or more fertility columns
#' @param stage a column name or position of the stage column in the stage-fate data frame. Defaults to "stage".
#' @param fate name of the fate column in the stage-fate data frame. Defaults to "fate". Individuals that die can be indicated with any character string not in the levels of stage, or missing.
#' @param fertility one or more names of fertility columns in the stage-fate data frame. By default, any column names matching stage class names are assumed to contain individual fertilities
#' @param sort a vector listing stage classes that correspond to the rows and columns of the desired projection matrix. Currently, names in this vector must match a level in the stage column. Also, this option should only be used if stages are not ordered, since the default is to sort by levels in the stage column.
#' @param add a vector listing row, column and value, used to add estimated transtions to the transition matrix (e.g., a transition from seed bank to seedling). May be repeated.
#' @param TF output separate transition (T) and fertility (F) matrices. Default is FALSE and outputs a single projection matrix A
#'
#' The state transition rates are estimated using transition frequency tables (see section 6.1.1, Caswell 2001), so this technique will most likely apply to demographic studies of plants or other sessile organisms where individuals are tagged and then consistently relocated in annual censuses. The fertility rates are calculated by averaging individuals fertilities by stage class; therefore, some care should be taken to correctly estimate individual fertilities based on the timing of the census.
#' Individual fertilities should be the total number of offspring at the end of the census interval. Therefore, fertilites should include offspring survival in a prebreeding censuses (and more than one offspring class may be present). In a postbreeding census, new offspring were born just before the census, so the fertility rate is just the number of offspring in this case.
#' @return The default output is a single projection matrix A. If the TF flag is true, then a list with 2 items where A=T+F
#' @export
#'
#' @source This is a modified version of \code{\link[popbio]{projection.matrix}}.
projection_matrix <- function (transitions, stage = NULL, fate = NULL, fertility = NULL,
          sort = NULL, add = NULL, TF = FALSE)
{
  if (missing(stage)) {
    stage <- "stage"
  }
  if (missing(fate)) {
    fate <- "fate"
  }
  nl <- as.list(1:ncol(transitions))
  names(nl) <- names(transitions)
  stage <- eval(substitute(stage), nl, parent.frame())
  fate <- eval(substitute(fate), nl, parent.frame())
  if (is.null(transitions[, stage])) {
    stop("No stage column matching ", stage)
  }
  if (is.null(transitions[, fate])) {
    stop("No fate column matching ", fate)
  }
  if (missing(sort)) {
    sort <- levels(transitions[, stage])
  }
  if (missing(fertility)) {
    fertility <- intersect(sort, names(transitions))
  }
  fertility <- eval(substitute(fertility), nl, parent.frame())
  tf <- table(transitions[, fate], transitions[, stage], useNA = "always")
  T_matrix <- try(prop.table(tf, 2)[sort, sort], silent = TRUE)
  if (inherits(T_matrix, "try-error")) {
    warning(paste("Error sorting matrix.\n  Make sure that levels in stage and fate columns\n  match stages listed in sort option above.\n Printing unsorted matrix instead!\n"),
            call. = FALSE)
    sort <- TRUE
    T_matrix <- prop.table(tf, 2)
  }
  T_matrix[is.nan(T_matrix)] <- 0
  if (length(add) > 0) {
    for (i in seq(1, length(add), 3)) {
      T_matrix[add[i + 0], add[i + 1]] <- as.numeric(add[i +
                                                           2])
    }
  }
  n <- length(fertility)
  F_matrix <- T_matrix * 0
  if (n == 0) {
    warning("Missing a fertility column with individual fertility rates\n",
            call. = FALSE)
  }
  else {
    for (i in 1:n) {
      fert <- tapply(transitions[, fertility[i]], transitions[,
                                                              stage], mean, na.rm = TRUE)[sort]
      F_matrix[i, ] <- fert
    }
  }
  F_matrix[is.na(F_matrix)] <- 0
  if (TF) {
    list(T = T_matrix, F = F_matrix)
  }
  else {
    T_matrix + F_matrix
  }
}
