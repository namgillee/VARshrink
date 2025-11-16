#' Restricted VAR
#'
#' Estimation of a VAR by imposing zero restrictions manually or by
#' significance.
#'
#' This is a modification of \code{vars::restrict()} for the class
#' "varshrinkest". Given an estimated VAR object of class "varest" or
#' "varshrinkest", a restricted VAR is obtained by choosing \code{method} "ser"
#' or "manual". Note: this function fits a restricted VAR using ordinary least
#' squares rather than a shrinkage method.
#'
#' @param x An object of class "varshrinkest"
#' @param ... Other arguments to vars::restrict()
#' @returns An object of class "varest", a VAR model fitted using ordinary
#' least squares. It contains an additional element \code{restrictions}, which
#' is a zero-one matrix indicating which variables were fixed to zero.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' restrict_sh(estim)
#' @seealso \code{\link[vars]{restrict}}
#' @export
restrict_sh <- function(x, ...) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  result <- vars::restrict(x, ...)
  return(result)
}
