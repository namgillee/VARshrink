print.varshrinksum <- function(varshsum) {
  # PRINT FUNCTIOn FOR THE CLASS "varshrinksum"
  #
  # See, summary.varshrinksum() for more about the class "varshrinksum"
  
  cat("Summary for VARshrink Estimation:\n")
  cat("=================================\n")
  
  if (is.null(varshsum$tsnames)) {
    mytsnames <- paste0("Y",1:d)
  } else {
    mytsnames <- varshsum$tsnames
  }
  cat("Endogenous variables: ")
  cat(mytsnames,"\n")
  cat("Deterministic variables: const")
  cat("Sample size: ", varshsum$obs)
  
  cat("log-likelihood:", varshsum$logLik, "\n")
  cat("AIC:", varshsum$AIC, "\n")
  cat("BIC:", varshsum$BIC, "\n")
  cat("Call:\n")
  print(varshsum$call)
  
  #---- FILL CODE HERE ----#
  cat("Estimated Coefficients:\n")
  # print cofficient matrix Psi (see, print.varshrinkest.R)
  
  cat("Standard Error of Coefficients:\n")
  # print SE of coefficient matrix
  
  cat("t value:\n")
  # print t value
  
  cat("Pr(>|t|):\n")
  # print p value
  
  cat("Signif. codes:  0 ¡®***¡¯ 0.001 ¡®**¡¯ 0.01 ¡®*¡¯ 0.05 ¡®.¡¯ 0.1 ¡® ¡¯ 1\n")
  # print signif. codes
  
  cat("\n")
  cat("Residual standard error: ")
  # Residual SE
  cat("on ", varshsum$df.residual, " degrees of freedom\n")
  
  cat("Multiple R-Squared: ")
  # multiple R-squared
  cat("       ")
  cat("Adjusted R-squared: ")
  # adjusted R-squared
  
  cat("\n")
  #cat("F-statistic: ")
  # skip F-statistic
  #cat(" on ")
  # numDF
  #cat(" and ")
  # denDF
  #cat(", ")
  #cat("p-value: ")
  # p-value
  #------------------------#
}