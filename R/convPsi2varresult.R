#' Convert format for VAR coefficients from Psi to varresult
#'
#' Convert a matrix of VAR coefficients estimated by a shrinkage method
#' into a list of "shrinklm" object. The class "shrinklm" inherits the class
#' "lm".
#'
#' Consider VAR(p) model: y_t = A_1 y_{t-1} + ... + A_p y_{t-p} + C d_t + e_t.
#' It can be written in the matrix form: Y = X %*% Psi + E,
#' where Psi is a concatenated M-by-K matrix, Psi = (A_1, ..., A_p, C)^T.
#'
#' The function converts Psi into a list of "shrinklm" objects which
#' inherited the class "lm". Consider the multiple linear regression form
#' of a VAR(p) model: y_j = X %*% psi_j + e_j, j=1,...,K.
#' Each "shrinklm" object contains the length-M vector psi_j as coefficients.
#'
#' Considering that each coefficient vector psi_j is estimated by a
#' shrinkage method, the effective number of parameters, k_eff, is
#' computed as: k_eff = Trace(X %*% (X^T %*% X + lambda0 * I)^{-1} %*% X').
#' Then, the degree of freedom of residuals is computed as: df.residual =
#' N - k_eff, where N is the number of rows of data matrices Y and X.
#'
#' @param Psi An M-by-K matrix of VAR coefficients
#' @param Y An N-by-K data matrix of dependent variables
#' @param X An N-by-M data matrix of regressors
#' @param lambda0 A rescaled shrinkage intensity parameter, based on which the
#' effective number of parameters is computed by Trace(X(X'X+lambda0*I)^{-1}X')
#' @param type Type of deterministic variables in the VAR estimation problem.
#' Either of "const", "trend", "both", or "none".
#' @param ybar,xbar NULL if Y and X are not centered. Mean vectors if Y and X
#' had been centered. If Y and X had been centered (ybar and xbar are not NULL)
#' and type is "const" or "both", then the coefficients for the constant term
#' is computed and concatenated to the coefficients.
#' @param Q_values Nonnegative weight vector of length N. Defaut is NULL.
#' Take weights on rows (samples) of Y and X by sqrt(Q).
#' @param callstr The call to VARshrink().
#' @return A list object with objects of class c("shrinklm", "lm").
#' Each "shrinklm" object has components: coefficients, residuals, fitted.values,
#' rank, df.residual, lambda0, call, terms, svd
#' @importFrom stats terms
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
