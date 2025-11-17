#' Causality Analysis
#'
#' Computes the test statistics for Granger- and Instantaneous causality for
#' a VAR(p).
#'
#' This function runs \code{vars::causality()} for an object of class
#' "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @param cause,vcov.,boot,boot.runs Other arguments for
#'   causality analysis; see help(causality) for details.
#' @returns A list of class attribute "htest"
#' with the following elements: \code{Granger}, \code{Instant}.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' causality_sh(estim, cause = "e")
#' @seealso \code{\link[vars]{causality}}
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
