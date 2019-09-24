#' Generate multivariate time series data using the given VAR model
#'
#' Generate a multivariate time series data set using the given VAR model.
#'
#' First, it creates (p+burnin+numT x K) data, then
#' it remove the first (p+burnin) vectors.
#' Finally, it returns (numT x K) data.
#'
#' @param numT Number of observed time points, T.
#' @param model A list object with Coef, Sigma, dof;
#' Coef is a list with A and c; A is a list object of K-by-K coefficient
#' matrices and c is a length-K vector.
#' Sigma is a K-by-K scale matrix and
#' dof is a degree of freedom for multivariate t-distribution for noise.
#' @param burnin Number of initial points which are not included in the final values.
#' @return A numT-by-K matrix
#' @examples
#' myCoef <- list(A = list(matrix(c(0.5, 0, 0, 0.5), 2, 2)), c = c(0.2, 0.7))
#' myModel <- list(Coef = myCoef, Sigma = diag(0.1^2, 2), dof = Inf)
#' simVARmodel(numT = 100, model = myModel, burnin = 10)
#' @export
simVARmodel <- function (numT, model, burnin = 0) {
  varCoef <- model$Coef
  noiseCov <- model$Sigma
  dof <- model$dof

  p <- length(varCoef$A)
  K <- dim(varCoef$A[[1]])[1] #K

  if (is.null(noiseCov)) {
    noiseCov <- diag(1, K) # Identity matrix
  } else {
    if ( length(noiseCov) == 1 ) {
      noiseCov <- diag(noiseCov, K) # diagonal matrix with equal diagonals
    }
  }

  #-- Generate All Noise Vectors --#
  if (is.infinite(dof)) {
    retTS <- MASS::mvrnorm(p + burnin + numT, rep(0, K), noiseCov)
  } else {
    retTS <- mvtnorm::rmvt(p + burnin + numT, df = dof, delta = rep(0, K), sigma = noiseCov)
  }

	#-- Generate time series --#

  #generate jth time series vector: burn_in period + T + p
	for (j in (p+1):(p+burnin+numT)) {
	    #const.vector
	    retTS[j, ] = retTS[j, ] + t(varCoef$c)

      #A(k) %*% y(t-k
	    for (k in 1:p) {
		      retTS[j, ] = retTS[j, ] + retTS[j - k, ] %*% t(varCoef$A[[k]])
	    }
	}

  colnames(retTS) <- paste("y", 1:K, sep = "")
	return( retTS[(p + burnin + 1):(p + burnin + numT), , drop = FALSE] )
}
