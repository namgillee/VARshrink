simVARmodel <- function (numT, model, burnin=0)

#simVARmodel <- function (numT, model, burnin=0)  
# Generate simulated multivariate time series data using the given VAR model
#
# First, it creates (p+burnin+numT x D) data, then 
# it remove the first (p+burning) vectors. 
# Finally, it returns (numT x D) data.
#
# Inputs  : 
#   numT : number of observed time points, T. 
#   model : an object of VARparam
#   burnin : number of initial points which are not included in the final values
#   
# Outputs : 
#   A matrix of size (numT x D)
#

{
  varCoef = model$Coef
  noiseCov = model$Sigma
  dof      = model$dof
  
  p = length(varCoef$A)
  d =  dim(varCoef$A[[1]])[1] #D
  
  if (is.null(noiseCov)) {
    noiseCov = diag(1,d) # Identity matrix
  } else { 
    if ( length(noiseCov)==1 ) {
      noiseCov = diag(noiseCov,d) # diagonal matrix with equal diagonals
    }
  }
  
  #-- Generate All Noise Vectors --#
  if (is.infinite(dof)) {
    require(MASS)
    retTS = mvrnorm(p+burnin+numT, rep(0,d), noiseCov)
  } else { 
    require(mvtnorm)
    retTS = rmvt(p+burnin+numT, df=dof, delta=rep(0,d), sigma=noiseCov)
  }
  
	#-- Generate time series --#

  #generate jth time series vector: burn_in period + T + p
	for (j in (p+1):(p+burnin+numT))
	{
	    #const.vector
	    retTS[j,] = retTS[j,] + t(varCoef$c)

      #A(k) %*% y(t-k
	    for (k in 1:p) 
	    {
		      retTS[j,] = retTS[j,] + retTS[j-k,] %*% t(varCoef$A[[k]])
	    }	                            
	}

	return( retTS[(p+burnin+1):(p+burnin+numT),] )
}