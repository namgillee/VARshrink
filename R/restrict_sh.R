#' Restrict function
#'
#' Warning: this function ignores shrinkage estimation, so it need to be
#' revised.
#'
#' @param x An object of class "varshrinkest"
#' @param ... Other arguments to vars::restrict()
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
