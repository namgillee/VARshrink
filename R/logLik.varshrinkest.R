#' Log-likelihood method for an object of class 'varshrinkest',
#' VAR parameters estimated by VARshrink()
#'
#' Returns the log-likelihood of an estimated VAR model
#' for class 'varshrinest'. Extend logLik.varest() to incorporate
#' 1) multivariate t-dribution for residuals,
#' 2) scale matrix Sigma provided by shrinkage methods, and
#' 3) effective number of parameters provided by shrinkage methods.
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University
# Acknowledgement: this code was contributed by Sung-Hoon Han & Dong-Han Lee
# @ Kangwon National University (2018.11.29.)
#' @export
logLik.varshrinkest <- function(object, ...) {

  obs <- object$obs
  df <- sum( obs - unlist(lapply(object$varresult, df.residual )) )
  K <- object$K
  resids <- resid(object) #
  Sigma <- if (is.null(object$Sigma)) {
    crossprod(resids) / obs
  } else {
    object$Sigma
  }
  dof <- ifelse(is.null(object$dof), Inf, object$dof)
  r <- if (is.infinite(dof)) {
    -(obs * K/2) * log(2 * pi) - (obs/2) * log(det(Sigma)) -
      (1/2) * sum( diag(resids %*% solve(Sigma, t(resids))) )
  } else {
    obs * lgamma((dof + K)/2) - obs * lgamma(dof / 2) -
      (obs * K/2) * log(dof * pi) - (obs/2) * log(det(Sigma)) -
      (dof + K)/2 * sum(log(1 +
                              diag(resids %*% solve(Sigma, t(resids))) / dof))
  }
  class(r) <- "logLik"
  attr(r, "df") <- df
  attr(r, "nobs") <- object$obs
  return(r)
}
