#' Normality, multivariate skewness and kurtosis test
#'
#' This function computes univariate and multivariate Jarque-Bera tests and
#' multivariate skewness and kurtosis tests for the residuals of a VAR(p)
#' or of a VECM in levels.
#' This is a modification of vars::normality.test() for
#' the class "varshrinkest".
#'
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @param multivariate.only If TRUE, only the multivariate statistics
#'   is computed.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' normality.test_sh(estim)
#' @seealso \code{\link[vars]{normality.test}}
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
