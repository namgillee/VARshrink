#' Restrict function
#'
#' Warning: this function ignores shrinkage estimation, so it need to be
#' revised.
#'
#' @param x An object of class "varshrinkest"
#' @param ... Other arguments to vars::restrict()
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
