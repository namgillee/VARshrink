lm_multiv_ridge <- function (Y, X, lambda = 0, do_scale = FALSE)
{
  #lm_multiv_ridge <- function (Y, X, lambda = 0, do_scale = FALSE)
  #  Multivariate ridge regression
  #
  # Inputs:
  #   Y        : Output data matrix Y having size N x K, or,
  #              a vector of length ND.
  #   X        : Input data matrix X having size N x M.
  #   do_scale : If true, center&scale X, and center Y.
  # Outputs:
  #   a list with the following attributes
  #      1) Psi - a list with as many Psi matrices as length(lambda)
  #      2) lambda - vector of lambda values
  #      3) GCV - vector of GCV values
  #

  p <- ncol(X)
  n <- nrow(X)
  if (!is.matrix(Y)) {
    dim(Y) <- c(n,length(Y)/n)
  }

  # Set lambda by a sequence of candidate values
  if (is.null(lambda)) {
    lambda = as.vector(c(1,5) %o% rep(10^c(-2:2)))
  }


  # Center and normalize X, center Y
  if (do_scale) {
    X <- scale(X, center = TRUE, scale = TRUE)
    Y <- scale(Y, center = TRUE, scale = FALSE)
  }

  # Compute SVD of X
  Xs <- svd(X)
  Rhs <- t(Xs$u) %*% Y        #r x K
  d <- Xs$d

  # Ridge regression: (X'X+nLI)^{-1} X'Y == V(d^2+nL)^{-1}d * Rhs
  k <- length(lambda)
  r <- length(d)
  Div <- d^2/n + rep(lambda, each=r*ncol(Y))   #(d^2/n + lambda)
  a <- rep(drop(d/n * Rhs), k)/Div             #(d^2/n + lambda)^{-1} * d/n * Rhs
  dim(a) <- c(r, ncol(Y)*k)
  coef <- Xs$v %*% a                           #p x K*k

  # GCV score
  Resid <- rep(Y,k) - X %*% coef
  dim(Resid) <- c(n*ncol(Y), k)                #n*K x k
  Divsmall <- d^2/n + rep(lambda, each=r)
  GCV <-  n * colSums(Resid^2) /  (n - colSums(matrix(d^2/n/Divsmall,r)))^2

  # Return list of M-by-K Psi matrices
  coef <- split(coef, rep(1:k, each=ncol(X)*ncol(Y)))
  coef <- lapply(coef, matrix, nrow = ncol(X), dimnames = list(colnames(X), colnames(Y)))

  res <- list(Psi = coef, lambda = lambda, GCV = GCV)
  res
}
