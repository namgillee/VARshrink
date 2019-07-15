#' @export
VARshrink  <- function(Y, p = 1, const_type = c('const', 'none', 'mean'),
                   method = c('ridge', 'ns', 'fbayes', 'sbayes','kcv'),
                   lambda = NULL, lambda_var = NULL, dof = Inf, ...)
# Shrinkage estimation of VAR parameters
#
#      y_t = sum_i A_i y_{t-i} + const. + eps_t
#
# Inputs:
#    Y:   T x D time series data
#    p:   lag order
#    const_type: 1) const - Estimate the const.
#                2) none  - Estimate without the const., i.e., c=0.
#                3) mean  - Use the sample mean, bar{y}, to center
#                     the time series data, and estimate the const.
#
#    method:     1) ridge - multivariate ridge regression
#                2) ns    - nonparametric shrinkage
#                3) fbayes- full Bayes MCMC shrinkage
#                4) sbayes- semi-parametric Bayes shrinkage
#                5) kcv   - k-fold cross validation
#
#    lambda, lambda_var:  shrinkage parameter value(s).
#                Use of this parameter is slightly different
#                for each method, that is, the same value does not
#                imply the same estimates.
#
#    dof:  degree of freedom of multivariate t distribution for noise.
#          Valid only for fbayes and sbayes.
#          If dof=Inf, it means multivariate normal distribution.
#
# Output:
#    an object of class 'varshrinkest' with components:
#       $varparam
#       $se.varparam
#       $p, $K, $obs, $totobs,
#       $const_type, $method,
#       $lambda, $lambda_var
#       $call
#
{
  cl <- match.call()
  rowY = nrow(Y) #T
  colY = ncol(Y) #D

  ######## Build Input-Output Data Matrix ########
  datX  = NULL
  datY  = NULL

  if (const_type == 'mean') {
    # Use the sample mean vector of Y
    ybar = colMeans(Y)
    Y = Y - rep(ybar, each=rowY)
  }

  datY = Y[-(1:p),]
  for (h in 1:p) {
    datX = cbind(datX, Y[(p+1-h):(rowY-h),])
  }
  N = nrow(datY)

  if (const_type == 'const') {
    # Append 1 to datX, because we will estimate
    # constant vector by shrinkage
    datX = cbind(1, datX)
  }


  ######## Run a Shrinkage Estimation Method ########

  estim = NULL

  #---------- (1) Multivariate Ridge ----------------
  if (method == 'ridge') {

    # Set lambda by a sequence
    if (is.null(lambda)) {
      lambda = as.vector(c(1,5) %o% rep(10^c(-2:2)))
    }

    # ridge estimates: $Coefs (list of all Psi_hat's), $lambda, $GCV
    resu_ridge  = lm_multiv_ridge(datY, datX, lambda, ...)

    # Select the minimum GCV
    id_min_gcv = which.min(resu_ridge$GCV)
    myPsi = resu_ridge$Coefs[[id_min_gcv]]

    # Residuals and covariance matrix for noise
    myResid = datY - datX %*% myPsi
    sing_val = svd(datX)$d #singular values
    kapp_val = sum( (sing_val^2)/(sing_val^2 + lambda[id_min_gcv]),
                    na.rm = TRUE) + identical(const_type,'mean') #Effective DOF
    mySigma = (t(myResid) %*% myResid)/max(1, N - kapp_val)

    # Convert Psi matrix into a list of $A and $c
    myCoef = convPsi2Coefs(myPsi)
    if (const_type == 'mean') {
      myCoef$c = convMean2Const(ybar, myCoef)
    }

    # Return the minimum GCV estimates
    estim$varparam = VARparam(Coef = myCoef, Sigma = mySigma, dof = Inf)
    estim$lambda = resu_ridge$lambda
    estim$lambda.estimated = (length(resu_ridge$lambda)>1)
    estim$df.residual = max(1, N - kapp_val)
    estim$GCV = resu_ridge$GCV

  }
  ##--------- (2) Nonparametric Shrinkage ----------
  if (method == 'ns') {
    require('corpcor') #Use cov.shrink() by Strimmer lab

    # Since NS method centers datX and datY,
    # We have to be careful to the case of const_type = 'const'.
    # Eventually, 'const' and 'mean' are identical:
    # If 'const'  -->  remove 1's from datX  -->  estimate c later by c(dybar, dxbar)
    # If 'mean'   -->         (no deletion)  -->  estimate c later by ybar
    # If 'none'   -->         (no deletion)  -->  (no estimation of c)

    if (const_type == 'const') {
      # Remove 1's from datX
      datX = datX[,-1]

      dybar = colMeans(datY)
      dxbar = colMeans(datX)

      datY = datY - rep(dybar, each=N)
      datX = datX - rep(dxbar, each=N)
    }

    # Estimate covariance matrix S_Z
    Z = cbind(datX, datY)
    if (is.null(lambda)) {
      if (is.null(lambda_var)){
        SZ = cov.shrink(Z, verbose = FALSE, ...)
      } else {
        SZ = cov.shrink(Z, lambda.var = lambda_var, verbose = FALSE, ...)
      }
    } else {
      if (is.null(lambda_var)){
        SZ = cov.shrink(Z, lambda = lambda, verbose = FALSE, ...)
      } else {
        SZ = cov.shrink(Z, lambda = lambda, lambda.var = lambda_var, verbose = FALSE, ...)
      }
    }

    # Compute the coefficient matrix Psi by solving linear system
    #if (attr(SZ,"lambda") > 1e-15) {
    #  myPsi = solve(SZ[1:(colY*p), 1:(colY*p)],
    #                  SZ[1:(colY*p), (colY*p+1):(colY*p+colY)])
    #} else {
    #  #S_X can be singular; Use eigen(S_X)
      eigSZ = eigen( SZ[1:(colY*p), 1:(colY*p)] )
      idr = eigSZ$values > 0
      myPsi = eigSZ$vectors[,idr] %*%
                ( 1/eigSZ$values[idr] * t(eigSZ$vectors[,idr]) ) %*%
                SZ[1:(colY*p), (colY*p+1):(colY*p+colY)]
    #}

    # Convert Psi matrix into a list of $A and $c
    myCoef = convPsi2Coefs(myPsi) #into list of A and c
    if (const_type == 'mean') {
      myCoef$c = convMean2Const(ybar, myCoef)
    }
    if (const_type == 'const') {
      # The colume of 1's has been removed from datX,
      # so the const_vector has to be estimated separately,
      # similarly to the case of const_type=='mean'.
      myCoef$c = matrix(dybar - t(myPsi) %*% dxbar,
                         nrow = colY, ncol = 1)
    }

    # Residuals and covariance matrix for noise
    myResid = datY - datX %*% myPsi
    sing_val = svd(datX)$d #singular values
    kapp_val = sum( (sing_val^2)/(sing_val^2 + attr(SZ,"lambda")*(N-1)/N),
                    na.rm = TRUE) + (const_type%in%c('mean','const')) #Effective DOF
    mySigma = (t(myResid) %*% myResid)/max(1, N - kapp_val)

    # Return the NS estimates
    estim$varparam = VARparam(Coef = myCoef, Sigma = mySigma, dof = Inf)
    estim$lambda = attr(SZ,"lambda")
    estim$lambda.estimated = attr(SZ,"lambda.estimated")
    estim$lambda_var = attr(SZ,"lambda.var")
    estim$lambda_var.estimated = attr(SZ,"lambda.var.estimated")
    estim$df.residual = max(1, N - kapp_val)
  }
  ##--------- (3) Full Bayesian with Noninformative Priors ----------
  if (method == 'fbayes') {

    # Estimate Psi coefficient matrix.
    # Arguments 'burnincycle' and 'mcmccycle' are included in '...'
    resu_fbayes = lm_full_Bayes_SR(datY, datX, dof = dof, ...)

    # Convert Psi matrix into a list of $A and $c
    myPsi = resu_fbayes$Psi
    myCoef = convPsi2Coefs(myPsi) #into list of A and c
    if (const_type == 'mean') {
      myCoef$c = convMean2Const(ybar, myCoef)
    }

    # Convert SE.Psi matrix into a list of $A and $c
    mySE.Psi = resu_fbayes$se.param$Psi
    mySE.Coef = convPsi2Coefs(mySE.Psi)
    if (const_type == 'mean') {
      tmpC01 =  convMean2Const( t(Y)%*%Y/(rowY-1), myCoef) #A*S
      mySE.Coef$c = diag( convMean2Const( t(tmpC01), myCoef) )#A*S'*A'=A*S*A'
      mySE.Coef$c = sqrt( mySE.Coef$c / (rowY) ) #sqrt( ASA'/n )
    }

    # Repeat for the LINEX estimator: Psi (optional)
    myPsi02 = resu_fbayes$LINEXparam$Psi
    myCoef02 = convPsi2Coefs(myPsi02)
    if (const_type == 'mean') {
      myCoef02$c = convMean2Const(ybar, myCoef02)
    }

    # Repeat for the LINEX estimator: SE.Psi (optional)
    mySE.Psi02 = resu_fbayes$se.LINEXparam$Psi
    mySE.Coefs02 = convPsi2Coefs(mySE.Psi02)
    if (const_type == 'mean') {
      tmpC02 =  convMean2Const( t(Y)%*%Y/(rowY-1), myCoef02) #A*S
      mySE.Coefs02$c = diag( convMean2Const( t(tmpC02), myCoef02) )#A*S'*A'=A*S*A'
      mySE.Coefs02$c = sqrt( mySE.Coefs02$c / (rowY) ) #sqrt( ASA'/n )
    }

    ## Return the fbayes estimates ##
    estim$varparam = VARparam(Coef = myCoef,
                              Sigma = resu_fbayes$Sigma,
                              dof = resu_fbayes$dof)
    estim$delta = resu_fbayes$delta
    estim$lambda = resu_fbayes$lambda
    estim$lambda.estimated = resu_fbayes$lambda.estimated
    estim$dof.estimated = resu_fbayes$dof.estimated

    estim$se.varparam = VARparam(Coef = mySE.Coef,
                                 Sigma = resu_fbayes$se.param$Sigma,
                                 dof = resu_fbayes$se.param$dof)
    estim$se.delta = resu_fbayes$se.param$delta
    estim$se.lambda = resu_fbayes$se.param$lambda

    estim$LINEXvarparam = VARparam(Coef = myCoef02) #An estimated coefficient minimizing LINEX loss
    estim$se.LINEXvarparam = VARparam(Coef = mySE.Coefs02) #The SE of the LINEX coefficient

    # df.residual
    # How to get kappa: From Eq. (11) for phi_hat(lambda), and defn. of hat matrix:
    #    Y = H*Y == X*Psi,
    #    ==> vec(Y) = y = (IxH)y == (IxX)psi
    #        and psi = (IxX'QX+lambda*SigxI)^{-1}(IxX'Q)y
    #    ==> Approximately, H == X(X'QX+lambda*s^2*I)^{-1}X'Q
    #        where  s^2 = tr(Sig)/D
    Q_values = resu_fbayes$q
    sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    s2_val = ifelse(colY>=2, mean(diag(resu_fbayes$Sigma), na.rm = TRUE), resu_fbayes$Sigma)
    kapp_val = sum( (sing_val^2)/(sing_val^2 + estim$lambda*s2_val/N),
                    na.rm = TRUE) + identical(const_type,'mean') #Effective DOF
    estim$df.residual = max(1, N - kapp_val)

  }
  ##--------- (4) Semi-parametric Bayesian with lambda by P-CV ----------
  if (method == 'sbayes') {

    resu_sbayes = lm_semi_Bayes_PCV(datY, datX, dof = dof,
                              lambda = lambda, lambda_var = lambda_var, ...)

    # Convert Psi matrix into a list of $A and $c
    myPsi = resu_sbayes$Psi
    myCoef = convPsi2Coefs(myPsi)
    if (const_type == 'mean') {
      myCoef$c = convMean2Const(ybar, myCoef)
    }

    ## Return the sbayes estimates
    estim$varparam = VARparam(Coef = myCoef,
                              Sigma = resu_sbayes$Sigma,
                              dof = resu_sbayes$dof)
    estim$lambda = resu_sbayes$lambda
    estim$lambda.estimated = resu_sbayes$lambda.estimated
    estim$lambda_var = resu_sbayes$lambda_var
    estim$lambda_var.estimated = resu_sbayes$lambda_var.estimated
    estim$dof.estimated = resu_sbayes$dof.estimated

    # df.residual
    Q_values = resu_sbayes$q
    sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    s2_val = ifelse(colY>=2, mean(diag(resu_sbayes$Sigma), na.rm = TRUE), resu_sbayes$Sigma)
    kapp_val = sum( (sing_val^2)/(sing_val^2 + estim$lambda/(1-estim$lambda)*s2_val*(N-1)/N),
                    na.rm = TRUE) + identical(const_type,'mean') #Effective DOF
    estim$df.residual = max(1, N - kapp_val)

  }
  ##--------- (5) Semi-parametric Bayesian with lambda by K-CV ----------
  if (method == 'kcv') {

    resu_kcv = lm_ShVAR_KCV(datY, datX, dof = dof,
                         lambda = lambda, lambda_var = lambda_var, ...)

    # Convert Psi matrix into a list of $A and $c
    myPsi = resu_kcv$Psi
    myCoef = convPsi2Coefs(myPsi) #into list of A's and c
    if (const_type == 'mean') {
      myCoef$c = convMean2Const(ybar, myCoef)
    }

    ## Return the kcv estimates
    estim$varparam = VARparam(Coef = myCoef,
                              Sigma = resu_kcv$Sigma,
                              dof = resu_kcv$dof)
    estim$lambda = resu_kcv$lambda
    estim$lambda.estimated = resu_kcv$lambda.estimated
    estim$lambda_var = resu_kcv$lambda_var
    estim$lambda_var.estimated = resu_kcv$lambda_var.estimated
    estim$dof.estimated = resu_kcv$dof.estimated

    # df.residual
    Q_values = resu_kcv$q
    sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    s2_val = ifelse(colY>=2, mean(diag(resu_kcv$Sigma), na.rm = TRUE), resu_kcv$Sigma)
    kapp_val = sum( (sing_val^2)/(sing_val^2 + estim$lambda/(1-estim$lambda)*s2_val*(N-1)/N),
                    na.rm = TRUE) + identical(const_type,'mean') #Effective DOF
    estim$df.residual = max(1, N - kapp_val)

  }
  ##---------------------------------------------------

  ## Check return value ##
  if (!with(estim, exists('varparam'))) {
    warning('VAR parameters were not estimated. Check the method.')
  }

  estim$residual = computeVARResidual(Y = Y, varparam = estim$varparam)
  estim$p        = p
  estim$K        = p*colY + 1 - identical(const_type,'none')  #ncol(datX) + identical(const_type,'mean')
  estim$obs      = N #rowY - p == nrow(datY) == N
  estim$totobs   = rowY
  estim$const_type = const_type
  estim$method   = method
  estim$call     = cl
  estim$tsnames  = colnames(Y)

  class(estim) <- c('varshrinkest')

  return(estim)
}
