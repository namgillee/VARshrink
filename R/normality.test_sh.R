#' Normality, multivariate skewness and kurtosis test for VAR shrinkage estimates
#'
#' An extension of vars::normallity.test() to the class "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @param multivariate.only If TRUE, only the multivariate statistics
#'   is computed.
#' @export
normality.test_sh <- function(x, multivariate.only = TRUE) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else if (inherits(x, "vec2var")) {
    class(x) <- "vec2var"
  } else {
    stop("\nPlease provide an object inheriting class 'varest' or 'vec2var'.\n")
  }
  result <- vars::normality.test(x, multivariate.only = multivariate.only)
  return(result)
}
