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
#' which is either 'const', or 'none'.
#' @param ybar,xbar mean vectors of Y and X which had been subtracted.
#' @param Q-values either NULL or a nonnegative weights of length N,
#' used for X'QX and X'QY in the hat matrix
#' @param callstr The call to VARshrink()
#' @return a list of c('shrinklm', 'lm') objects with components:
#' coefficients, residuals, fitted.values, rank, df.residual, lambda0,
#' call, terms, svd

convPsi2varresult <- function(Psi, Y, X, lambda0,
                              type = c('const', 'none'),
                              ybar = NULL, xbar = NULL,
                              Q_values = NULL,
                              callstr = "")
{
  N = nrow(Y)
  K = ncol(Y)
  p = nrow(Psi) %/% K

  #=========================================================================#
  # Compute the fitted values and the residuals: myFitted, myResid (N-by-K)
  #
  # myFitted = X %*% Psi
  # myResid = Y - myFitted
  #=========================================================================#
  myFitted = X %*% Psi
  myResid = Y - myFitted

  dof_to_adjust = 0 #number of parameters to adjust
  ##-------- REMOVE THE OPTION type=='mean' --------##
  # if (identical(tolower(type), 'mean')) {
  #   #----------------------------------------------
  #   # If type=='mean', then
  #   #  1) Compute and append the const vector to Psi, and
  #   #  2) Add ybar to the fitted values.
  #   # Assume that:
  #   #     ybar == colMeans(y),
  #   #     X == {x_t' - ybar'},
  #   #     Y == {y_t' - ybar'}.
  #   # Compute:
  #   #   const' = ybar' - ybar' %*% Psi
  #   # MODIFY myFitted:
  #   #   yhat_t' = const' + x_t' %*% Psi
  #   #           = ybar' + (x_t' - ybar') %*% Psi
  #   #           = ybar' + X %*% Psi .
  #   # DO NOT NEED TO MODIFY myResid:
  #   #  r_t' = y_t' - yhat_t'
  #   #       = y_t' - ybar' - X %*% Psi
  #   #       = Y - X %*% Psi
  #   #-------------------------------------------
  #   if (!is.null(ybar)) {
  #     if (length(ybar) == K) {
  #       const = as.vector(ybar) - as.vector(t(rep(ybar, p)) %*% Psi)
  #       Psi = rbind(Psi, const = const) #1) Append the const vector to Psi
  #       dof_to_adjust = 1
  #
  #       myFitted = myFitted + rep(ybar, each = N)  #2) Add ybar to the fitted values
  #
  #     } else {
  #       warning("Length of ybar is incorrect")
  #     }
  #
  #   } else {
  #     # Assume that ybar is a zero vector
  #     Psi = rbind(Psi, const = rep(0, K)) #1) Append the const vector to Psi
  #     dof_to_adjust = 1
  #
  #     #2) Add ybar to the fitted values: skip
  #   }
  # }
  ##------------------------------------------------##
  if (identical(tolower(type), 'const') &&
      (nrow(Psi) %% K == 0) &&
      !is.null(ybar) && !is.null(xbar)) {
    #---------------------------------------------
    # If type=='const' but Psi does not include the
    # const vector in its last row, then suppose
    # that method=='ns' and
    #  1) Compute and append the const vector to Psi, and
    #  2) Add ybar to the fitted values.
    # Assume that
    #      ybar == colMeans(datY),
    #      xbar == colMeans(datX).
    # Compute    const = ybar - xbar %*% Psi;
    # because    yhat_t' = const' + x_t' %*% Psi
    # and the 'ns'  equation    Y = X %*% Psi
    # with Y=={yhat_t'-ybar'}, X=={x_t'-xbar'}.
    # Add ybar to the fitted values because
    #           yhat_t' = const' + x_t' %*% Psi
    #                   = const' + X %*% Psi + xbar' %*% Psi
    #                   = ybar + X %*% Psi
    #---------------------------------------------

    const = as.vector(ybar) - as.vector(t(xbar) %*% myPsi)
    Psi = rbind(Psi, const = const) #1) Append the const vector to Psi
    dof_to_adjust = 1

    myFitted = myFitted + rep(ybar, each = N) #2) Add ybar to the fitted values
  }

  #=========================================================================#
  # Compute the effective number of parameters: mykapp
  #   based on Trace(X(X'X+lambda0*I)^{-1}X')
  #=========================================================================#
  if (is.null(Q_values)) {
    s = svd(X)
    sing_val = s$d
  } else {
    s = svd(sqrt(Q_values)*X)
    sing_val = s$d
  }
  idnonzero <- (abs(sing_val) >= 1e-14)
  mykapp = sum((sing_val[idnonzero]^2)/(sing_val[idnonzero]^2 + lambda0), na.rm = TRUE) +
    dof_to_adjust
  #=========================================================================#
  # Return value
  #=========================================================================#
  varresult = vector('list', K)
  names(varresult) = colnames(Y)
  obj = list(coefficients = Psi[,1],
             residuals = myResid[,1],
             fitted.values = myFitted[,1],
             rank = sum(idnonzero) + dof_to_adjust,
             df.residual = max(1, N - mykapp),
             lambda0 = lambda0,
             call = callstr,
             terms = terms(y ~ ., data = X),
             svd = s # for class 'shrinklm', we replace qr=qr(X) with svd=svd(X)
  )
  class(obj) <- c('shrinklm', 'lm')

  for(i in 1:K){
    obj$coefficients = Psi[,i]
    obj$residuals = myResid[,i]
    obj$fitted.values = myFitted[,i]
    varresult[[i]] <- obj
  }

  return(varresult)

}
