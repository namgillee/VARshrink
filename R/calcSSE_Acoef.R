#' Sum of squared error (SSE) between two VAR coefficients
#'
#' VAR coefficients in form of a list of A matrices.
#'
#' @param Acoef1,Acoef2 list of coefficient matrices
#' @return SSE value

#' @export
calcSSE_Acoef <- function(Acoef1, Acoef2) {

	p1 = length(Acoef1)
	p2 = length(Acoef2)
	if (p1 != p2)  {
	  warning("Input Acoef1 and Acoef2 have different orders..")
	  return(-1)
	}
	sseValue <- 0
	for (i in 1:p1) {
	  sseValue = sseValue + sum((Acoef1[[i]] - Acoef2[[i]])^2)
	}
	return(sseValue)
}

