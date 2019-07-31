#' Summary method for an object of class 'varshrinkest',
#' VAR parameters estimated by VARshrink()
#'
#' Extend summary.varest() to class 'varshrinest' to incorporate
#' adapted methods for new classes:
#' summary.shrinklm(), logLik.varshrinkest(), roots.varshrinkest().
#'
#' Code is modified to avoid call to data matrices ($y, $datamat)
#' and to use effective numbers of parameters of shrinkage estimates.
#'
#' Output includes the scale matrix, Sigma, and degree-of-freedom, dof,
#' for multivariate t-distribution for residuals.
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University

summary.varshrinkest <- function (object, equations = NULL, ...) {
  ynames <- names(object$varresult)
  obs <- object$obs
  if (is.null(equations)) {
    ysubnames <- ynames
  }
  else {
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
  roots <- roots.varshrinkest(object)
  result <- list(names = ysubnames, varresult = eqest, covres = covres,
                 corres = corres, logLik = logLik, obs = obs, roots = roots,
                 type = object$type, call = object$call,
                 Sigma = Sigma, dof = dof)
  class(result) <- c("varshrinksum", "varsum")
  return(result)
}
