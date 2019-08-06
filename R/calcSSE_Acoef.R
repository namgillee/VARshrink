#' Sum of squared errors (SSE) between coefficients of two VARs
#'
#' Compute sum of squared errors of coefficients of lagged endogenous
#' variables (Acoef) of two VAR models.
#'
#' Consider VAR(p) model: y_t = A_1 y_{t-1} + ... + A_p y_{t-p} + C d_t + e_t.
#' The SSE of two VAR(p) models is expressed as
#' sum_{k=1}^p sum_{i=1}^K sum_{j=1}^K ( (A_k)_{ij} - (A_k')_{ij} )^2.
#'
#' @param Acoef1,Acoef2 Each one is a list object with K-by-K coefficient
#' matrices of lagged endogenous variables. See help(Acoef.varshrinkest), or,
#' help(Acoef).
#' @return SSE value.
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

