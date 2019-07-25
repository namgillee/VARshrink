#' Convert Psi to varresult
#'
#' Consider the VAR process, Y = X %*% Psi + E, in the matrix form.
#' This function returns a list of lm objects for the multiple
#' regression problems corresponding to each column of Y:
#' y_j = X %*% psi_j + e_j.
#'
#' @param Psi M-by-K matrix Psi
#' @param Y N-by-K matrix
#' @param X N-by-M matrix
#' @param lambda a shrinkage intensity parameter, which will be used
#' for computing the effective number of parameters
#' @param scale_lambda a factor to be multiplied to lambda
#' @param type Either 'const', 'mean', or 'none', indicating the
#' constant type in the VAR shrinkage problem.
#' @param ybar global mean vector if type=='mean'; local mean if type=='const'
#' @param xbar local mean vector for datX if type=='const'
#' @param callstr The call to VARshrink()
#' @return a list of 'lm' objects with components: coefficients,
#' residuals, fitted.values, rank, df.residual, call

convPsi2varresult <- function(Psi, Y, X, lambda, scale_lambda = 1,
                              type = c('const', 'mean', 'none'),
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
  if (identical(tolower(type), 'mean')) {
    #-------------------------------------------------------
    # If type=='mean', then assume that:
    #     ybar == colMeans(y),
    #     X == {x_t' - ybar'},
    #     Y == {y_t' - ybar'}.
    # Fitted value for the VAR equation: const' = ybar' - ybar' %*% Psi
    #  yhat_t' = const' + x_t' %*% Psi
    #          = ybar' + (x_t' - ybar') %*% Psi
    #          = ybar' + X %*% Psi .
    #     (NEED TO MODIFY myFitted BY ADDING ybar)
    # Residual is:
    #  r_t' = y_t' - yhat_t'
    #       = y_t' - ybar' - X %*% Psi
    #       = Y - X %*% Psi
    #     (DO NOT NEED TO MODIFY myResid)
    #-------------------------------------------------------
    if (!is.null(ybar)) {
      if (length(ybar) == K) {
        const = as.vector(ybar) - t(rep(ybar, p)) %*% Psi
        Psi = rbind(Psi, const) #2) Append the const vector to Psi
        dof_to_adjust = 1

        myFitted = myFitted + rep(ybar, each = N)  #1) Add ybar to the fitted values

      } else {
        warning("Length of ybar is incorrect")
      }

    } else {
      # Assume that ybar is a zero vector
      Psi = rbind(Psi, const = rep(0, K)) #2) Append the const vector to Psi
      dof_to_adjust = 1

      #1) Add ybar to the fitted values: skip
    }
  }
  if (identical(tolower(type), 'const') &&
      (nrow(Psi) %% K == 0) &&
      !is.null(ybar) && !is.null(xbar)) {
    #-------------------------------------------------------
    # In the case that type=='const' but Psi does not include
    # the const vector in one of its rows, assume that
    #      ybar == colMeans(datY),
    #      xbar == colMeans(datX),
    # we consider that the case is for method=='ns', so
    #   1) Compute
    #           const = ybar - xbar %*% Psi;
    #      because the VAR equantion is
    #           yhat_t' = const' + x_t' %*% Psi
    #      and the 'ns' estimation equatin is Y = X %*% Psi
    #      with Y=={yhat_t'-ybar'}, X=={x_t'-xbar'}
    #   2) Psi = rbind(Psi, const)
    #   3) Add ybar to the fitted values because
    #           yhat_t' = const' + x_t' %*% Psi
    #                   = const' + X %*% Psi + xbar' %*% Psi
    #                   = ybar + X %*% Psi
    #-------------------------------------------------------

    const = as.vector(ybar) - t(xbar) %*% myPsi
    Psi = rbind(Psi, const) #2) Append the const vector to Psi
    dof_to_adjust = 1

    myFitted = myFitted + rep(ybar, each = N) #1) Add ybar to the fitted values
  }

  #=========================================================================#
  # Compute the effective number of parameters: mykapp
  #   based on Trace(X(X'X+lambda*I)^{-1}X')
  #=========================================================================#
  if (is.null(Q_values)) {
    sing_val = svd(X)$d #singular values of X
  } else {
    sing_val = svd(sqrt(Q_values)*X)$d #singular values
  }
  mykapp = sum((sing_val^2)/(sing_val^2 + lambda * scale_lambda), na.rm = TRUE) +
    dof_to_adjust

  # Return value
  varresult = vector('list', K)
  names(varresult) = colnames(Y)
  obj = list(coefficients = Psi[,1],
             residuals = myResid[,1],
             fitted.values = myFitted[,1],
             rank = sum(abs(sing_val) > 1e-14),
             df.residual = max(1, N - mykapp),
             call = callstr
  )
  class(obj) <- 'lm'

  for(i in 1:K){
    obj$coefficients = Psi[,i]
    obj$residuals = myResid[,i]
    obj$fitted.values = myFitted[,i]
    varresult[[i]] <- obj
  }

  return(varresult)

}
