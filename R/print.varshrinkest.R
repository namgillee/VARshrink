print.varshrinkest <- function(varshest) {
  # PRINT FUNCTION FOR THE CLASS "varshrinkest"
  
  cat("VARshrink Estimation Results:\n")
  cat("=============================\n")
  cat("Call:\n")
  print(varshest$call)
  
  #### TS variable names ####
  d = nrow(varshest$varparam$Coef$A[[1]])
  p = length(varshest$varparam$Coef$A)
  if (is.null(varshest$tsnames)) {
    mytsnames <- paste0("Y",1:d)
  } else {
    mytsnames <- varshest$tsnames
  }
  
  #### VAR Coefficients ####
  cat("\n")
  cat("Estimated coefficients:\n")
  Psi = matrix(unlist(varshest$varparam$Coef$A), d)
  rownames(Psi) = mytsnames
  colnames(Psi) = paste0(mytsnames, ".lag", rep(1:p,each=d) )
  if (varshest$const_type != "none") {
    Psi = cbind(Psi, varshest$varparam$Coef$c)
    colnames(Psi)[d*p+1] = "const"
  } 
  print(round(Psi, 6))
  
  #### Other Parameters ####
  cat("\n")
  cat("$Sigma for noise:\n")
  Sig <- varshest$varparam$Sigma
  rownames(Sig) <- colnames(Sig) <- mytsnames
  print(round(Sig,6))
  
  cat("\n")
  cat("$dof for noise:", varshest$varparam$dof,
      "(estimated:", ifelse(is.null(e<-varshest$dof.estimated),FALSE,e),
      ")","\n")
  cat("$lambda:", varshest$lambda,
      "(estimated:", ifelse(is.null(e<-varshest$lambda.estimated),FALSE,e),
      ")","\n")
  #### Optional Parameters ####
  if (!is.null(varshest$GCV)) { ## For ridge
    cat("$GCV:", round(varshest$GCV,5),"\n")
  }
  if (!is.null(varshest$lambda_var)) { ## For NS and sbayes
    cat("$lambda_var:", varshest$lambda_var,
        "(estimated:", ifelse(is.null(e<-varshest$lambda_var.estimated),FALSE,e),
        ")","\n")
  }
} 