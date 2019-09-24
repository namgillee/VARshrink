#' Test for serially correlated errors for VAR shrinkage estimate
#'
#' An extension of vars::serial.test() to the class "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @param lags.pt,lags.bg,type Other arguments for vars::serial.test().
#'   see help(serial.test) for details.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' serial.test_sh(estim)
#' @seealso \code{\link[vars]{serial.test}}
#' @export
serial.test_sh <- function(x, lags.pt = 16, lags.bg = 5, type =
                             c("PT.asymptotic",
                               "PT.adjusted", "BG", "ES")) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else if (inherits(x, "vec2var")) {
    class(x) <- "vec2var"
  } else {
    stop("\nPlease provide an object inheriting class 'varest' or class 'vec2var'.\n")
  }
  result <- vars::serial.test(x, lags.pt = lags.pt,
                              lags.bg = lags.bg, type = type)
  return(result)
}
