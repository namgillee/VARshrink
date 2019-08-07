#' Restrict function
#'
#' Warning: this function ignores shrinkage estimation, so it need to be
#' revised.
#'
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
