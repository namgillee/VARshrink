#' @export
VARparam      <- function(Coef, Sigma = NULL, dof = Inf) {
  # The set of VAR model parameters
  # The VAR Model equation is
  #      y(t) = A[[1]]*y(t-1) + ... + A[[p]]*y(t-p) + c + e,
  # with e having a multivarate t distribution with
  # the covariance matrix 'Sigma' and degree-of-freedom 'dof'

  resu = list(Coef = Coef, Sigma = Sigma, dof = dof)
  attr(resu, 'class') <- 'VARparam'
  resu
}
