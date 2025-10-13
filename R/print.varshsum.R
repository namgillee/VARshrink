#' Print method for class "varshsum"
#'
#' Print method for an object obtained by summary.varshrinkest().
#'
#' This function extends print.varsum() for VAR models estimated by
#' shrinkage methods.
#' The output includes scale matrix Sigma and degree of freedom dof
#' for multivariate t-distribution for residuals.
#'
# Last modified: 2019.7.30. Namgil Lee @ Kangwon National University
#' @param x An object of class "varshsum"
#' @param digits,signif.stars,... Other arguments for print(),
#' printCoefmat(), format() method
#' @importFrom stats printCoefmat pf cov2cor
#' @export
print.varshsum <- function(x, digits = max(3, getOption("digits") - 3),
                           signif.stars = getOption("show.signif.stars"),
                           ...) {
  dim <- length(x$names)
  text1 <- "\nVAR Shrinkage Estimation Results:\n"
  cat(text1)
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")
  cat(paste("Endogenous variables:",
            paste(colnames(x$covres),
                  collapse = ", "), "\n", collapse = " "))
  cat(paste("Deterministic variables:", paste(x$type, collapse = ", "),
            "\n", collapse = " "))
  cat(paste("Sample size:", x$obs, "\n"))
  cat(paste("Log Likelihood:", round(x$logLik, 3), "\n"))
  cat("Roots of the characteristic polynomial:\n")
  cat(formatC(x$roots, digits = digits))
  cat("\nCall:\n")
  print(x$call)
  cat("\n\n")
  for (i in 1:dim) {
    result <- x$varresult[[x$names[i]]]
    text1 <- paste("Estimation results for equation ", x$names[i],
                   ":", sep = "")
    cat(text1, "\n")
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    text2 <- paste(x$names[i], " = ", paste(rownames(result$coef),
                                            collapse = " + "), sep = "")
    cat(text2, "\n\n")
    printCoefmat(result$coef, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    cat("\n")
    cat("\nResidual standard error:",
        format(signif(result$sigma, digits)), "on", result$df[2L],
        "degrees of freedom\n")
    if (!is.null(result$fstatistic)) {
      cat("Multiple R-Squared:", formatC(result$r.squared,
                                         digits = digits))
      cat(",\tAdjusted R-squared:", formatC(result$adj.r.squared,
                                            digits = digits),
          "\nF-statistic:", formatC(result$fstatistic[1], digits = digits),
          "on", result$fstatistic[2],
          "and", result$fstatistic[3], "DF,  p-value:",
          format.pval(pf(result$fstatistic[1L], result$fstatistic[2L],
                         result$fstatistic[3L], lower.tail = FALSE),
                      digits = digits), "\n")
    }
    cat("\n\n")
  }
  cat("\nScale matrix, Sigma, of multivariate t distribution for noise:\n")
  print(x$Sigma, digits = digits, ...)
  cat("\nDegrees of freedom of multivariate t distribution for noise:\n")
  print(x$dof, digits = digits, ...)
  cat("\nCorrelation matrix of Sigma:\n")
  print(cov2cor(x$Sigma), digits = digits, ...)
  cat("\n\n")
  invisible(x)
}
