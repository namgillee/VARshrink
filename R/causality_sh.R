#' Causality Analysis for class "varshrinkest"
#'
#' An extension of vars::causality() to the class "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @param cause,vcov.,boot,boot.runs Other arguments for
#'   causality analysis; see help(causality) for details from
#'   the package vars documentation.
#' @export
causality_sh <- function(x, cause = NULL, vcov. = NULL,
                         boot = FALSE, boot.runs = 100) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  x$datamat <- as.data.frame(x$datamat)
  result <- vars::causality(x, cause = cause, vcov. = vcov., boot = boot,
                            boot.runs = boot.runs)
  return(result)
}
