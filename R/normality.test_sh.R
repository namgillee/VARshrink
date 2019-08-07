#' @export
normality.test_sh <- function(x, ...) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else if (inherits(x, "vec2var")) {
    class(x) <- "vec2var"
  } else {
    stop("\nPlease provide an object inheriting class 'varest' or 'vec2var'.\n")
  }
  result <- vars::normality.test(x, ...)
  return(result)
}
