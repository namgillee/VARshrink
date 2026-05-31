#' Summary method for an object of class "varshrinkest",
#' VAR parameters estimated by VARshrink()
#'
#' Extend \code{summary.varest()} to the class "varshrinkest" to incorporate
#' adapted methods for the new classes:
#' \code{summary.shrinklm()}, \code{logLik.varshrinkest()}.
#'
#' The code has been modified to eliminate direct calls to the data matrices
#' (\code{$y}, \code{$datamat}) and to use the effective numbers of parameters
#' obtained from the shrinkage estimates. The output additionally includes the
#' covariance matrix \code{Sigma} and the degrees-of-freedom \code{dof}
#' of the multivariate t-distribution for the residuals.
#'
#' @param object An object of class "varshrinkest", usually
#' a result of call to "VARshrink()".
#' @param equations Subset of names of endogenous time series variables
#' to summarize.
#' @param ... Currently not used.
#' @importFrom stats resid cov df.residual cor
#' @returns Returns a list with class attribute "varshsum" and "varsum" which
#' contains the following elements: names, logLik, obs, roots, type, call,
#' varresult, covres, corres, Sigma, dof.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' summary(estim)
#' @export
summary.varshrinkest <- function(object, equations = NULL, ...) {
  ynames <- names(object$varresult)
  obs <- object$obs
  if (is.null(equations)) {
    ysubnames <- ynames
  } else {
    ysubnames <- as.character(equations)
    if (!(all(ysubnames %in% ynames))) {
      warning("\nInvalid variable name(s) supplied, using first variable.\n")
      ysubnames <- ynames[1]
    }
  }
  eqest <- lapply(object$varresult[ysubnames], summary)
  resids <- resid(object)
  Sigma <- if (is.null(object$Sigma)) {
    crossprod(resids) / obs
  } else {
    object$Sigma
  }
  dof <- ifelse(is.null(object$dof), Inf, object$dof)
  covres <- cov(resids) * (obs - 1) / min(sapply(object$varresult, df.residual))
  corres <- cor(resids)
  logLik <- as.numeric(logLik(object))
  roots <- roots(object)
  result <- list(names = ysubnames, varresult = eqest, covres = covres,
                 corres = corres, logLik = logLik, obs = obs, roots = roots,
                 type = object$type, call = object$call,
                 Sigma = Sigma, dof = dof)
  class(result) <- c("varshsum", "varsum")
  return(result)
}
