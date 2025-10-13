#' ARCH-LM test
#'
#' Performs univariate and multivariate ARCH-LM tests for a VAR.
#' This is a modification of vars::arch.test() for the class "varshrinkest".
#'
#' @param x An object of class "varshrinkest" obtained by VARshrink()
#' @param lags.single An integer of the lag order used for
#'   univariate ARCH statistics.
#' @param lags.multi An integer of the lag order used for
#'   multivariate ARCH statistic.
#' @param multivariate.only If TRUE, only the multivariate statistic
#'   is computed.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' arch.test_sh(estim)
#' @seealso \code{\link[vars]{arch.test}}
#' @export
arch.test_sh <- function(x, lags.single = 16, lags.multi = 5,
                         multivariate.only = TRUE) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else if (inherits(x, "vec2var")) {
    class(x) <- "vec2var"
  } else {
    stop(paste0("\nPlease provide an object inheriting class 'varest' or class",
                "'vec2var'.\n"))
  }
  result <- vars::arch.test(x, lags.single = lags.single,
                            lags.multi = lags.multi,
                            multivariate.only = multivariate.only)
  return(result)
}
