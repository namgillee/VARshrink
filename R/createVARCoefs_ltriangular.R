createVARCoefs_ltriangular <- function(p=1, d=5, diag_val=1/p, num_nonzero=0, const_vector = NULL, range_min=0.2, range_max=1/p)
#createVARCoefs_ltriangular <- function(p=1, d=5, diag_val=1/p, num_nonzero=0, const_vector = NULL, range_min=0.1, range_max=1/p)
# Create coefficient matrices A1, ..., Ap and a constant vector c for a VARmodel where
#   diag(A_i) = range_max
# and other nonzero coefficients are selected from
# lower-triangular parts of A_i, whose values
# are drawn from a uniform distribution Unif(\pm [range_min, range_max])
# 
# Inputs:
#  p : lag order
#  d : number of time series variables $D$
#  diag_val : diagonal values of A1,...,Ap
#  num_nonzero : number of nonzero entries on the lower-triangular parts of  A1,..,Ap
#  const_vector : constant vector $c$ in the VAR model
#  range_min, range_max : Each nonzero entry of coefficient matrices
#      is drawn uniformly from the interval
#      [-range_max, -range_min] U [range_min, range_max]
# 
# Outputs:
#  A list of $A and $c, which is suitable for VARparam$Coef. 
#  $A is a list of matrices A1, ..., Ap, 
#  and $c is a constant vector. 
#  (See, VARparam$Coef)
#


{
	#-- Generate coefficients ---#
  varCoef = vector('list',2)
  names(varCoef) = c('A','c')
  
  coef = matrix(0,d,d*p)
	for (j in 1:p) 
  {
    diag(coef[,(1+(j-1)*d):(j*d)]) = diag_val
	}
	
  #-- Draw nonzero coefficients ---%
  if (floor(num_nonzero) > 0) 
  {
    #INDEX FOR LOWER-TRIANGULAR PARTS OF EACH A_i,i=1,..,p
    idxzero <- NULL
    lower_indices <- which(lower.tri(matrix(0,d,d)))
    for (j in 1:p) {
          idxzero <- c(idxzero, lower_indices + (j - 1)*d*d)
    }    
    lenzero <- length(idxzero)
    
    idx_nonzero <- sample(lenzero,  min(lenzero,floor(num_nonzero)) ) #selec location
  	val_nonzero = runif(length(idx_nonzero), min=-1, max=1) #values
  	val_nonzero[val_nonzero >= 0] = (range_max - range_min) * val_nonzero[val_nonzero >= 0] + range_min
  	val_nonzero[val_nonzero < 0] = (range_max - range_min) * val_nonzero[val_nonzero < 0] - range_min
    
  	coef[idxzero[idx_nonzero]] = val_nonzero
  }
  
  #-- Return value --#
	for (j in 1:p) {
		varCoef$A[[j]] = coef[,(1+(j-1)*d):(j*d)]
	}
	
	if (is.null(const_vector)) {
		varCoef$c = matrix(0, d, 1)
	}	else 	{
		varCoef$c = const_vector
	}

	return(varCoef)
}