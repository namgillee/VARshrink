#' Structural stability of a VAR(p)
#'
#' Computes empirical fluctuation processes for VAR estimates.
#' Utilizes \code{strucchange::efp()} for the VAR estimates of each time series
#' variable.
#'
#' A variant of \code{vars::stability()} for an object of class "varshrinkest".
#' @param x Object of class "varshrinkest"
#' @param type,h,dynamic,rescale,... Other arguments to
#' \code{strucchange::efp()}
#' @returns A list with class attribute "varstabil" which contains the following
#' elements: \code{stability}, \code{endog}, \code{K}. The \code{stability} is
#' a list of \code{strucchange::efp()} outputs.
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
#' stabil <- stability_sh(estim)
#' plot(stabil)
#' @seealso \code{\link[vars]{stability}}
#' @export
stability_sh <- function(x, type = c("OLS-CUSUM", "Rec-CUSUM", "Rec-MOSUM",
                                     "OLS-MOSUM", "RE", "ME", "Score-CUSUM",
                                     "Score-MOSUM", "fluctuation"),
                         h = 0.15, dynamic = FALSE, rescale = TRUE, ...) {
  type <- match.arg(type)
  K <- x$K
  stability <- list()
  endog <- colnames(x$datamat)[1:K]
  for (i in 1:K) {
    formula <- formula(x$varresult[[i]])
    data <- cbind(x$datamat[, endog[i]], x$datamat[, -c(1:K)])
    colnames(data)[1] <- "y"
    data <- as.data.frame(data)
    stability[[endog[i]]] <-
      strucchange::efp(formula = formula, data = data,
                       type = type, h = h, dynamic = dynamic,
                       rescale = rescale)
  }
  result <- list(stability = stability, names = endog, K = K)
  class(result) <- "varstabil"
  return(result)
}
