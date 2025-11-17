#' Multivariate Ridge Regression
#'
#' Estimate regression coefficients by using ridge regression.
#'
#' Consider the multivariate regression:
#' \deqn{\mathbf{Y} = \mathbf{X \Psi} + \mathbf{e}.}
#' \eqn{\mathbf{\Psi}} is a M-by-K matrix of regression coefficients.
#' The ridge regression estimate for the coefficients is
#' \deqn{\mathbf{\Psi} = (\mathbf{X'X} + \lambda \mathbf{I})^{-1}
#' \mathbf{X'Y}.}
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param lambda Numeric vector of lambda values
#' @param do_scale If true, X is centered and scaled, and Y is centered.
#' @return A list object with the components: 1) Psi - A list of
#' estimated \eqn{\mathbf{\Psi}} matrices, 2) lambda - A vector of
#' \eqn{\lambda} values, 3) GCV - A vector of GCV values
#' @references G. H. Golub, M. Heath, G. Wahba (1979).
#' Generalized cross-validation as a method for choosing a good
#' ridge parameter. Technometrics 21(2), 215-223. doi: 10.2307/1268518

lm_multiv_ridge <- function(Y, X, lambda = 0, do_scale = FALSE) {

  n <- nrow(X)
  if (!is.matrix(Y)) {
    dim(Y) <- c(n, length(Y) / n)
  }

  # Set lambda by a sequence of candidate values
  if (is.null(lambda)) {
    lambda <- as.vector(c(1, 5) %o% rep(10 ^ c(-4:1)))
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
  Div <- d ^ 2 / n + rep(lambda, each = r * ncol(Y))
  a <- rep(drop(d / n * Rhs), k) / Div
  dim(a) <- c(r, ncol(Y) * k)
  coef <- Xs$v %*% a                           #p x K*k

  # GCV score
  Resid <- rep(Y, k) - X %*% coef
  dim(Resid) <- c(n * ncol(Y), k)                #n*K x k
  Divsmall <- d ^ 2 / n + rep(lambda, each = r)
  GCV <-  n * colSums(Resid ^ 2) /
    (n - colSums(matrix(d ^ 2 / n / Divsmall, r))) ^ 2

  # Return list of M-by-K Psi matrices
  coef <- split(coef, rep(1:k, each = ncol(X) * ncol(Y)))
  coef <- lapply(coef, matrix, nrow = ncol(X),
                 dimnames = list(colnames(X), colnames(Y)))

  res <- list(Psi = coef, lambda = lambda, GCV = GCV)
  res
}
