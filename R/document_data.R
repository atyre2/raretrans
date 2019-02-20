#' Add roxygen2 documentation template for a dataset
#'
#' @param object data set to document
#'
#' @keywords internal
document_data <- function(object){
  x1 <- "#' ~~ data name/kind ... \n#'"
  x1b <- "#' ~~ A concise (1-5 lines) description of the dataset. ~~\n#'"
  x2 <- "#' \\itemize{"
  x3 <- paste(paste0("#'    \\item ", names(object), "..."), collapse = "\n")
  x4 <- "#'  }\n#'"
  x5 <- "#' @docType data"
  x6 <- "#' @keywords ..."
  x7 <- paste("#' @name", deparse(substitute(object)))
  x8 <- paste0("#' @usage data(", deparse(substitute(object)), ")")
  x9 <- paste("#' @format A data frame with", NROW(object), "rows and",
              NCOL(object), "variables"
  )
  x = paste(c(x1, x1b, x2, x3, x4, x5, x6, x7, x8, x9), collapse = "\n")
  writeLines(x, con = file.path("R", paste0(deparse(substitute(object)), '.R')))
}
