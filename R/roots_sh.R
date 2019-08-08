#' Eigenvalues of the companion coefficient matrix of
#' a VAR(p)-process
#'
#' This is a variant of vars::roots() for an object of class 'varshrinkest',
#' VAR parameters estimated by VARshrink()
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University
#' @param x An object of class "varshrinkest"
#' @param modulus TRUE for modulus of the roots.
#' @export
roots_sh <- function(x, modulus = TRUE) {

  if (!inherits(x, "varest")) {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  K <- x$K
  p <- x$p
  A <- unlist(Acoef_sh(x))
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  roots <- eigen(companion)$values
  if (modulus)
    roots <- Mod(roots)
  return(roots)
}
