#' BQ function for class "varshrinkest"
#'
#' An extension of vars::BQ() to the class "varshrinkest".
#' @param x An object of class "varshrinkest" obtained by VARshrink().
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' BQ_sh(estim)
#' @seealso \code{\link[vars]{BQ}}
#' @export
BQ_sh <- function(x) {
  if (inherits(x, "varest")) {
    class(x) <- "varest"
  } else {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  Amats <- Acoef_sh(x)
  P <- x$p
  Ident <- diag(x$K)
  mat1 <- matrix(0, x$K, x$K)
  mat2 <- mat1
  for (i in 1:P) {
    mat1 <- mat1 - Amats[[i]]
    mat2 <- mat2 - t(Amats[[i]])
  }
  mat1 <- Ident + mat1
  mat2 <- Ident + mat2
  df <- summary(x$varresult[[1]])$df[2]
  SigmaU <- crossprod(resid(x))/df
  eval <- solve(mat1) %*% SigmaU %*% solve(mat2)
  lrim <- t(chol(eval))
  colnames(lrim) <- colnames(x$y)
  rownames(lrim) <- colnames(lrim)
  cim <- mat1 %*% lrim
  colnames(cim) <- colnames(lrim)
  rownames(cim) <- colnames(lrim)
  result <- list(A = Ident, Ase = NULL, B = cim, Bse = NULL,
                 LRIM = lrim, Sigma.U = SigmaU * 100, LR = NULL, opt = NULL,
                 start = NULL, type = "Blanchard-Quah", var = x,
                 call = match.call())
  class(result) <- "svarest"
  return(result)
}
