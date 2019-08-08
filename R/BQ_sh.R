#' BQ function for class "varshrinkest"
#'
#' An extension of vars::BQ() to the class "varshrinkest".
#' @export
BQ_sh <- function(x) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  result <- vars::BQ(x)
  return(result)
}
