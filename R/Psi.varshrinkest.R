Psi.varshrinkest <- function (x, nstep = 10, ...) {
  if (!(inherits(x, "varest"))) {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  nstep <- abs(as.integer(nstep))
  Phi <- Phi(x, nstep = nstep)
  Psi <- array(0, dim = dim(Phi))
  dfr <- df.residual(x$varresult[[1]])
  sigma.u <- crossprod(resid(x)) / dfr
  P <- t(chol(sigma.u))
  dim3 <- dim(Phi)[3]
  for (i in 1:dim3) {
    Psi[, , i] <- Phi[, , i] %*% P
  }
  return(Psi)
}
