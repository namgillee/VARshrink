#' Convert Psi to varresult
#'
#' Consider the VAR process, Y = X %*% Psi + E, in the matrix form.
#' This function returns a list of 'shrinklm' objects, which is
#' a slight modification to the 'lm', for the multiple
#' regression problems corresponding to each column of Y:
#' y_j = X %*% psi_j + e_j.
#'
#' Depending on the 'type', the fitted values (y_j) and the estimated
#' coefficients (psi_j) are adjusted automatically.
#'
#' @param Psi M-by-K matrix Psi
#' @param Y N-by-K matrix
#' @param X N-by-M matrix
#' @param lambda0 a rescaled shrinkage intensity parameter
#' to be used for computing the effective number of parameters
#' @param type the type of the constant vector in the VAR shrinkage problem,
#' which is one of "const", "trend", "both", or "none".
#' @param ybar,xbar mean vectors of Y and X which had been subtracted.
#' @param Q-values either NULL or a nonnegative weights of length N,
#' used for X'QX and X'QY in the hat matrix
#' @param callstr The call to VARshrink()
#' @return a list of c('shrinklm', 'lm') objects with components:
#' coefficients, residuals, fitted.values, rank, df.residual, lambda0,
#' call, terms, svd

convPsi2varresult <- function(Psi, Y, X, lambda0,
                              type = c("const", "trend", "both", "none"),
                              ybar = NULL, xbar = NULL,
                              Q_values = NULL,
                              callstr = "") {

  N <- nrow(Y)
  K <- ncol(Y)

  #### Compute the fitted values and the residuals ####

  my_fitted <- X %*% Psi
  my_resid <- Y - my_fitted

  nparam_to_adjust <- 0 #number of parameters to adjust

  if ((identical(type, "const") || identical(type, "both")) &&
      !is.null(ybar) && !is.null(xbar)) {
    # If both ybar and xbar are supplied as arguments, then
    # suppose that:   X = (x_t - xbar), Y = (y_t - ybar), Y ~ X %*% Psi,
    # and
    # 1) estimate the constant vector by const = ybar - Psi' %*% xbar,
    # 2) append const to Psi, and
    # 3) add ybar to the fitted values. (Because:
    #    yhat_t' = const' + x_t' %*% Psi = ybar + X %*% Psi.)

    const <- as.vector(ybar) - as.vector(t(xbar) %*% Psi)
    Psi <- rbind(Psi, const = const)
    nparam_to_adjust <- nparam_to_adjust + 1

    my_fitted <- my_fitted + rep(ybar, each = N)
  }

  #### Compute the effective number of parameters: mykapp ####
  ####   based on Trace(X(X'X+lambda0*I)^{-1}X')          ####

  if (is.null(Q_values)) {
    s <- svd(X)
    sing_val <- s$d
  } else {
    s <- svd(sqrt(Q_values) * X)
    sing_val <- s$d
  }
  idnonzero <- (abs(sing_val) >= 1e-14)
  mykapp <- sum( (sing_val[idnonzero]^2) /
                   (sing_val[idnonzero]^2 + lambda0),
                 na.rm = TRUE) + nparam_to_adjust

  #### Return value ####

  varresult <- vector("list", K)
  names(varresult) <- colnames(Y)
  obj <- list(coefficients = Psi[, 1],
              residuals = my_resid[, 1],
              fitted.values = my_fitted[, 1],
              rank = sum(idnonzero) + nparam_to_adjust,
              df.residual = max(1, N - mykapp),
              lambda0 = lambda0,
              call = callstr,
              terms = terms(y ~ ., data = X),
              svd = s # for class 'shrinklm', replace qr=qr(X) with svd=svd(X)
  )
  class(obj) <- c("shrinklm", "lm")

  for (i in 1:K){
    obj$coefficients <- Psi[, i]
    obj$residuals <- my_resid[, i]
    obj$fitted.values <- my_fitted[, i]
    varresult[[i]] <- obj
  }

  return(varresult)

}
