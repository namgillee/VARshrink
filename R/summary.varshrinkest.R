logLik.varshrinkest <- function(varshest) {
  # Compute log-likelihood of the estimated VAR parameters
  #
  # The input object varshest is an object of the class "varshrinkest" 
  # obtained from VARshrink.R
  #
  # 2018.11.29. Sung-Hoon Han & Dong-Han Lee & Namgil Lee @ Kangwon National University


  # Res is residuals defined by: dataY - datX %*% Psi
  # Sig and dof are covariance matrix and DOF 
  #for multivariate t-distribution for noise 
  Res <- varshest$residual    # Residuals 
  Sig <- varshest$varparam$Sigma  # Estimated covariance matrix
  dof <- varshest$varparam$dof  # Estimated DOF
  
  if (is.null(dof)) {
    stop("In logLik():dof is NULL")
  }
  
  if (is.infinite(dof)) {
    # Use multivariate normal distribution
    #c.f.: retTS = mvrnorm(p+burnin+numT, rep(0,d), noiseCov)
    #---------- FILL THIS CODE ---------------#
    require(mvtnorm)
    return(sum(dmvnorm(Res, mean=rep(0,ncol(Res)), Sig, log=TRUE)))
    #-----------------------------------------#
    
  } else { 
    # Use multivariate t distribution
    
    #c.f.: retTS = rmvt(p+burnin+numT, df=dof, delta=rep(0,d), sigma=noiseCov)
    #---------- FILL THIS CODE ---------------#
    require(mvtnorm)
    return(sum(dmvt(Res, delta=rep(0,ncol(Res)), sigma=Sig, df=dof, log=TRUE)))
    #-----------------------------------------#
  }
}

AIC.varshrinkest <- function(varshest, k=2) {
  # AIC
  #
  # 2018.11.23. Namgil Lee @ Kangwon National Univeristy
  
  N <- varshest$obs
  
  kapp <- N - varshest$df.residual  ## Effective DOF
  
  return( -2 * logLik(varshest) - k * kapp )
}

BIC.varshrinkest <- function(varshest) {
  # BIC
  #
  # 2018.11.23. Namgil Lee @ Kangwon National Univeristy
  
  return( AIC.varshrinkest(varshest, k=log(varshest$obs)) )
}

summary.varshrinkest <- function(varshest) {
  # SUMMARY FUNCTION FOR THE CLASS "varshrinkest"
  # Returns an object of class "varshrinksum"
  #
  # See, VARshrink.R for more about the class "varshrinkest"
  #
  # 2018.11.23. Namgil Lee @ Kangwon National University
  
  ret <- NULL
  
  ret$varparam <- varshest$varparam #estimated parameters
  
  ret$logLik <- logLik(varshest) #log-likelihood
  
  ret$AIC <- AIC(varshest)
  
  ret$BIC <- BIC(varshest)
  
  ret$df.residual <- varshest$df.residual
  
  ret$obs <- varshest$obs
  
  ret$totobs <- varshest$totobs
  
  #ret$roots
  
  ret$const_type <- varshest$const_type
  
  ret$call <- varshest$call
  
  ret$tsnames  <- varshest$tsnames
  
  class(ret) <- "varshrinksum"
  return(ret)
  
}
