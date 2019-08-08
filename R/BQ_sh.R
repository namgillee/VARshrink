#' BQ function for class "varshrinkest"
#'
#' An extension of vars::BQ() to the class "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
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
