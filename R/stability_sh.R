#' Stability function
#'
#' Warning: this function has not been tested for small sample sizes yet.
#'
#' @export
stability_sh <- function (x, type = c("OLS-CUSUM", "Rec-CUSUM", "Rec-MOSUM",
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
