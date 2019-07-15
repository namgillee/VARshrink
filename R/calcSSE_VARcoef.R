#' @export
calcSSE_VARcoef <- function(model1, model2, include_const_vector = FALSE)

#Calculate sum of squared error (SSE) to measure
#the difference between coefficients of two VARs
#
# Inputs:
#     model1, model2 :  objects of VARparam
# Outputs:
#     MSE value

{
  d = dim(model1$Coef$A[[1]])[1]
	p1 = length(model1$Coef$A)
	p2 = length(model2$Coef$A)


	if (p1 != p2)  {
    warning("Input model1 and model2 have different orders..")
    return(-1)
	}


	#sseValue <- 0
	# (1) constant vector
	if (include_const_vector) {
	  #num_elem <- d*d*p1+d
	  sseValue <- sum((model1$Coef$c - model2$Coef$c)^2)
	} else {
	  #num_elem = d*d*p1
	  sseValue <- 0
	}
	# (2) coefficients
	for (i in 1:p1) {
	  sseValue = sseValue + sum((model1$Coef$A[[i]] - model2$Coef$A[[i]])^2)
	}

	return(sseValue)
}

