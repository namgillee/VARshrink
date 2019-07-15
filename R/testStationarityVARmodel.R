#' @export
testStationarityVARmodel <- function(model)

#testStationarityVARmodel <- function(model)
# Returns TRUE if the given VAR model is covariance stationary,
# based on Hamilton (1994), proposition 10.1.
#
# Inputs:
#  model : an object of VARparam
#
# Outputs:
#  TRUE or FALSE
#
#------
# The VAR Model is covariance stationary if |x| < 1 for all
# values of x satisfying
#
#      det(F-xI) = | I x^p - A1 x^{p-1} - ... - Ap | = 0.
#
# In other words, all eigenvalues of the marix F is |x|<1.
#------

{
  varCoef = model$Coef

	p = length(varCoef$A)
	d = dim(varCoef$A[[1]])[1]

  #-- matrix F -----#
	matrixF = varCoef[[1]] # the matrix F
  if (p>=2)
	{
            for (k in 2:p)
                matrixF <- cbind(matrixF, varCoef$A[[k]])

            matrixB <- cbind( diag(1,(p-1)*d), matrix(0,(p-1)*d,d) )

            matrixF <- rbind(matrixF, matrixB)
  }

  #-- result --#

  if ( max(abs(eigen(matrixF)$values)) >= 1 )
	    return (FALSE)
	else
	    return (TRUE)

}
