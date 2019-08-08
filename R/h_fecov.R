#' @importFrom stats df.residual resid
h_fecov <- function (x, n.ahead) {
  n.par <- sapply(x$varresult, df.residual)
  sigma.u <- crossprod(resid(x))/n.par
  Sigma.yh <- array(NA, dim = c(x$K, x$K, n.ahead))
  Sigma.yh[, , 1] <- sigma.u
  Phi <- Phi(x, nstep = n.ahead)
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp <- matrix(0, nrow = x$K, ncol = x$K)
      for (j in 2:i) {
        temp <- temp + Phi[, , j] %*% sigma.u %*% t(Phi[, , j])
      }
      Sigma.yh[, , i] <- temp + Sigma.yh[, , 1]
    }
  }
  return(Sigma.yh)
}
