#' Eigenvalues of the companion coefficient matrix of
#' a VAR(p)-process
#'
#' This is a variant of vars::roots() for an object of class 'varshrinkest',
#' VAR parameters estimated by VARshrink()
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University
#' @export
roots_sh <- function(x, ...) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  result <- vars::roots(x, ...)
  return(result)
}
