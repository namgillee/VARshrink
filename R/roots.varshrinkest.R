#' Eigenvalues of the companion coefficient matrix of
#' a VAR(p)-process
#'
#' This is a variant of vars::roots() for an object of class 'varshrinkest',
#' VAR parameters estimated by VARshrink()
#'
#' Extend the roots() to class 'varshrinest' to incorporate
#' adapted methods: Acoef.varshrinkest(), Bcoef.varshrink
#'
#' Code is modified to avoid call to data matrices ($y, $datamat)
#' and to use effective numbers of parameters of shrinkage estimates.
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University

roots.varshrinkest <- function(x, modulus = TRUE) {
  if (!inherits(x, "varest")) {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  K <- x$K
  p <- x$p
  A <- unlist(Acoef.varshrinkest(x))
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
