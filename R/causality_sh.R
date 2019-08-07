#' @export
causality_sh <- function(x, ...) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  x$datamat <- as.data.frame(x$datamat)
  result <- vars::causality(x, ...)
  return(result)
}
