#' Coefficient matrices (A coefficients) for an object of class 'varshrinkest'
#'
#' A variant of vars::Acoef(). Code is modified to avoid call to data matrices
#' ($y, $datamat) and to use effective numbers of parameters of shrinkage
#' estimates.
#' @export
Acoef.varshrinkest <- function (x) {

  if (!inherits(x, "varest")) {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  K <- x$K
  p <- x$p
  A <- Bcoef.varshrinkest(x)[, 1:(K * p)]
  As <- list()
  start <- seq(1, p * K, K)
  end <- seq(K, p * K, K)
  for (i in 1:p) {
    As[[i]] <- matrix(A[, start[i]:end[i]], nrow = K, ncol = K)
    rownames(As[[i]]) <- rownames(A)
    colnames(As[[i]]) <- colnames(A[, start[i]:end[i]])
  }
  return(As)
}
