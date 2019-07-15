#' @export
computeVARResidual <- function(Y, varparam) {

  # Compute residuals from multivariate time series data and
  # estimated VAR parameters: Roughly,
  #   Res = datY - sum_i(datYi %*% t(Ai)) - c
  #
  # Input
  #   Y          : TxD time series data
  #   varparam   : object of class "VARparam"
  #
  # Output
  #   Res        : (T-p)xD matrix of residuals

  p <- length(varparam$Coef$A)
  #d <- nrow(varparam$Coef$A[[1]])
  lenT <- nrow(Y)

  # Compute residuals
  Res <- Y[-(1:p),]
  for (i in 1:p) {
    Res <- Res - Y[(p+1-i):(lenT-i),] %*% t(varparam$Coef$A[[i]])
  }
  if (!is.null(varparam$Coef$c))
    Res <- Res - rep(varparam$Coef$c, each=lenT-p)

  return(Res)
}
