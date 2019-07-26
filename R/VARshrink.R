#' Shrinkage estimation of VAR parameters
#'
#' @param y T-by-K time series data
#' @param p lag order
#' @param type   1) const - Estimate the const.
#                2) mean  - Use the sample mean, bar{y}, to center
#                     the time series data, and estimate the const.
#                3) none  - Estimate without the const., i.e., c=0.
#
#' @param method 1) ridge - multivariate ridge regression
#'                2) ns    - nonparametric shrinkage
#'                3) fbayes- full Bayes MCMC shrinkage
#'                4) sbayes- semi-parametric Bayes shrinkage
#'                5) kcv   - k-fold cross validation
#'
#' @param lambda,lambda_var  shrinkage parameter value(s).
#'                Use of this parameter is slightly different
#'                for each method, that is, the same value does not
#'                imply the same estimates.
#'
#' @param dof  degree of freedom of multivariate t distribution for noise.
#'              Valid only for fbayes and sbayes.
#'             If dof=Inf, it means multivariate normal distribution.
#' @return an object of class c('varshrinkest','varest') including the following fields:
#' varresult, datamat, y, type, p, K, obs, totobs, restrictions, method,
#' lambda, lambda_var, call
#'
#' @importFrom corpcor cov.shrink
#' @import vars
#' @export
VARshrink  <- function(y, p = 1, type = c('const', 'mean', 'none'),
                   method = c('ridge', 'ns', 'fbayes', 'sbayes','kcv'),
                   lambda = NULL, lambda_var = NULL, dof = Inf, ...)

{
  y <- as.matrix(y)
  totobs = nrow(y)   #total number of observations
  K = ncol(y)        #dimension of output response
  N = totobs - p     #sample size
  M = K * p + ifelse(identical(tolower(type), "const"), 1, 0) #dimension of input covariates

  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if (K < 2)
    stop("The matrix 'y' should contain at least two variables.\n")
  if (is.null(tsnames <- colnames(y)))
    tsnames <- paste("y", 1:K, sep = "")

  if (totobs <= (p+1) )
    stop("Number of total observations must be > p+1\n")

  ######## Build Data Matrices: datX, datY ########

  #  i) ybar is computed when datY, datX are centered by colMeans(y)
  # ii) dybar,dxbar are computed when datY, datX are centered by colMeans(datY), colMeans(datX)
  ybar = NULL
  dybar = dxbar = NULL
  if (identical(tolower(type), 'mean')) {
    #********** Subtract the mean from y **********
    ybar = colMeans(y)
    y = y - rep(ybar, each = totobs)
  }

  datY = y[(p+1):totobs, ] #N-by-K
  colnames(datY) <- tsnames
  datX = matrix(1, N, M) #N-by-M
  for (h in 1:p) {
    # Note that, if type=='const', then 1's are at the last column, since M>(K*p)
    datX[, (1+(h-1)*K):(h*K)] = y[(p+1-h):(totobs-h), ]
  }
  if(identical(tolower(type), "const")) {
    #datX has M=K*p+1 columns
    colnames(datX) <- c(paste(rep(tsnames, times = p) , ".l", rep(1:p, each = K), sep = ""), "const")
  } else {
    #datX has M=K*p columns
    colnames(datX) <- paste(rep(tsnames, times = p) , ".l", rep(1:p, each = K), sep = "")
  }

  #### Run a Shrinkage Estimation Method: estim ####

  estim = NULL

  #---------- (1) Multivariate Ridge ----------------
  if (method == 'ridge') {

    # Compute the coefficient matrix: myPsi (M-by-K)
    # resu_ridge is a list of $Psi, $lambda, $GCV
    resu_ridge  = lm_multiv_ridge(datY, datX, lambda, ...)
    id_min_gcv = which.min(resu_ridge$GCV) # Select the minimum GCV
    myPsi = resu_ridge$Psi[[id_min_gcv]]

    # Update the return value
    estim$varresult <- convPsi2varresult(Psi = myPsi, Y = datY, X = datX,
                                         lambda = lambda[id_min_gcv], scale_lambda = 1,
                                         type = type, ybar = ybar,
                                         callstr = cl
    )
    estim$lambda = resu_ridge$lambda
    estim$lambda.estimated = as.logical(length(resu_ridge$lambda) > 1)
    estim$GCV = resu_ridge$GCV

  }
  ##--------- (2) Nonparametric Shrinkage ----------
  if (method == 'ns') {

    # In 'ns', datX and datY are centered separately, by dxbar and dybar.
    # If type=='const', remove 1's from datX, and estimate Psi, and
    #  estimate const later by const' = dybar' - dxbar' %*% Psi.
    # If type=='mean' or 'none', no need to remove 1's from datX.

    if (identical(tolower(type), 'const')) {
      datX = datX[,-M]  # Remove 1's from datX

      dybar = colMeans(datY)
      dxbar = colMeans(datX)

      datY = datY - rep(dybar, each = N)
      datX = datX - rep(dxbar, each = N)
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
    eigSZ = eigen(SZ[1:(K*p), 1:(K*p)]) #Use EVD to avoid singular matrix problem
    idr = eigSZ$values > 0
    myPsi = eigSZ$vectors[,idr] %*%
            ( 1/eigSZ$values[idr] * t(eigSZ$vectors[,idr]) ) %*%
            SZ[1:(K*p), (K*p+1):(K*p+K)]
    rownames(myPsi) <- colnames(datX); colnames(myPsi) <- colnames(datY)

    # Update the return value
    estim$varresult <- convPsi2varresult(Psi = myPsi, Y = datY, X = datX,
                                         lambda = attr(SZ,"lambda"), scale_lambda = (N-1)/N,
                                         type = type, ybar = dybar, xbar = dxbar,
                                         callstr = cl
    )
    estim$lambda = attr(SZ,"lambda")
    estim$lambda.estimated = attr(SZ,"lambda.estimated")
    estim$lambda_var = attr(SZ,"lambda.var")
    estim$lambda_var.estimated = attr(SZ,"lambda.var.estimated")

  }
  ##--------- (3) Full Bayesian with Noninformative Priors ----------
  if (method == 'fbayes') {

    # Compute the coefficient matrix: myPsi
    # Arguments 'burnincycle' and 'mcmccycle' are included in '...'
    resu_fbayes = lm_full_Bayes_SR(datY, datX, dof = dof, ...)

    # Update the return value
    estim$lambda = resu_fbayes$lambda
    estim$lambda.estimated = resu_fbayes$lambda.estimated
    estim$delta = resu_fbayes$delta

    estim$Sigma = resu_fbayes$Sigma
    estim$dof = resu_fbayes$dof
    estim$dof.estimated = resu_fbayes$dof.estimated

    myPsi = resu_fbayes$Psi
    rownames(myPsi) <- colnames(datX); colnames(myPsi) <- colnames(datY)
    sigbar = ifelse(K>=2, mean(diag(resu_fbayes$Sigma), na.rm = TRUE), resu_fbayes$Sigma)
    estim$varresult <- convPsi2varresult(Psi = myPsi, Y = datY, X = datX,
                                         lambda = resu_fbayes$lambda, scale_lambda = sigbar,
                                         type = type, ybar = ybar,
                                         Q_values = resu_fbayes$q,
                                         callstr = cl
    )

    # ## overwrite df.residual ##
    # Q_values = resu_fbayes$q
    # sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    # sigbar = ifelse(K>=2, mean(diag(resu_fbayes$Sigma), na.rm = TRUE), resu_fbayes$Sigma)
    # mykapp = sum((sing_val^2)/(sing_val^2 + resu_fbayes$lambda*sigbar),
    #              na.rm = TRUE) + as.vector(identical(type,'mean'))
    # for (i in 1:K) {
    #   estim$varresult[[i]]$df.residual = max(1, N - mykapp)
    # }

    # # Convert SE.Psi matrix into a list of $A and $c
    # mySE.Psi = resu_fbayes$se.param$Psi
    # mySE.Coef = convPsi2Coefs(mySE.Psi)
    # if (type == 'mean') {
    #   tmpC01 =  convMean2const( t(y)%*%y/(totobs-1), myCoef) #A*S
    #   mySE.Coef$c = diag( convMean2const( t(tmpC01), myCoef) )#A*S'*A'=A*S*A'
    #   mySE.Coef$c = sqrt( mySE.Coef$c / (totobs) ) #sqrt( ASA'/n )
    # }
    #
    # # Repeat for the LINEX estimator: Psi (optional)
    # myPsi02 = resu_fbayes$LINEXparam$Psi
    # myCoef02 = convPsi2Coefs(myPsi02)
    # if (type == 'mean') {
    #   myCoef02$c = convMean2const(ybar, myCoef02)
    # }
    #
    # # Repeat for the LINEX estimator: SE.Psi (optional)
    # mySE.Psi02 = resu_fbayes$se.LINEXparam$Psi
    # mySE.Coefs02 = convPsi2Coefs(mySE.Psi02)
    # if (type == 'mean') {
    #   tmpC02 =  convMean2const( t(y)%*%y/(totobs-1), myCoef02) #A*S
    #   mySE.Coefs02$c = diag( convMean2const( t(tmpC02), myCoef02) )#A*S'*A'=A*S*A'
    #   mySE.Coefs02$c = sqrt( mySE.Coefs02$c / (totobs) ) #sqrt( ASA'/n )
    # }
    #
    # estim$se.varparam = VARparam(Coef = mySE.Coef,
    #                              Sigma = resu_fbayes$se.param$Sigma,
    #                              dof = resu_fbayes$se.param$dof)
    # estim$se.lambda = resu_fbayes$se.param$lambda
    # estim$se.delta = resu_fbayes$se.param$delta
    #
    # estim$LINEXvarparam = VARparam(Coef = myCoef02) #An estimated coefficient minimizing LINEX loss
    # estim$se.LINEXvarparam = VARparam(Coef = mySE.Coefs02) #The SE of the LINEX coefficient

    #### How to get kappa ####
    # From Eq. (11) for psi_hat(lambda), Y can be separated as:
    #    psi := psi_hat(lambda)
    #        = (I(x)X'QX + lambda*Sig(x)I)^{-1} %*% (I(x)X'Q) %*% y
    # From defn. of hat matrix:
    #    Y_hat = H %*% Y == X %*% Psi
    # After vectorizing:
    #    vec(Y_hat) = (I(x)H) %*% y == (I(x)X) %*% psi
    # Combining above two equations:
    #    (I(x)H) %*% y == (I(x)X) %*% (I(x)X'QX + lambda*Sig(x)I)^{-1} %*% (I(x)X'Q) %*% y
    #    ...
    #    (I(x)H) == (I(x)X) %*% (I(x)X'QX + lambda*Sig(x)I)^{-1} %*% (I(x)X'Q)
    # Approximate Sig by a scalar:
    #    Sig == sigbar := Tr(Sig)/K
    # then,
    #    H == X %*% (X'QX + lambda * Sigbar * I)^{-1} %*% X'Q
    #    where  sigbar := Tr(Sig)/K

  }
  ##--------- (4) Semi-parametric Bayesian with lambda by P-CV ----------
  if (method == 'sbayes') {

    resu_sbayes = lm_semi_Bayes_PCV(datY, datX, dof = dof,
                              lambda = lambda, lambda_var = lambda_var, ...)

    # Update the return value
    estim$lambda = resu_sbayes$lambda
    estim$lambda.estimated = resu_sbayes$lambda.estimated
    estim$lambda_var = resu_sbayes$lambda_var
    estim$lambda_var.estimated = resu_sbayes$lambda_var.estimated

    estim$Sigma = resu_sbayes$Sigma
    estim$dof = resu_sbayes$dof
    estim$dof.estimated = resu_sbayes$dof.estimated

    myPsi = resu_sbayes$Psi
    rownames(myPsi) <- colnames(datX); colnames(myPsi) <- colnames(datY)
    sigbar = ifelse(K>=2, mean(diag(resu_sbayes$Sigma), na.rm = TRUE), resu_sbayes$Sigma)
    estim$varresult <- convPsi2varresult(Psi = myPsi, Y = datY, X = datX,
                                         lambda = resu_sbayes$lambda, scale_lambda = sigbar,
                                         type = type, ybar = ybar,
                                         Q_values = resu_sbayes$q,
                                         callstr = cl
    )

    ## overwrite df.residual ##
    # Q_values = resu_sbayes$q
    # sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    # sigbar = ifelse(K>=2, mean(diag(resu_sbayes$Sigma), na.rm = TRUE), resu_sbayes$Sigma)
    # mykapp = sum((sing_val^2)/(sing_val^2 + resu_sbayes$lambda/(1-resu_sbayes$lambda)*sigbar*(N-1)/N),
    #              na.rm = TRUE) + as.vector(identical(type,'mean'))
    # for (i in 1:K) {
    #   estim$varresult[[i]]$df.residual = max(1, N - mykapp)
    # }

  }
  ##--------- (5) Semi-parametric Bayesian with lambda by K-CV ----------
  if (method == 'kcv') {

    resu_kcv = lm_ShVAR_KCV(datY, datX, dof = dof,
                         lambda = lambda, lambda_var = lambda_var, ...)

    # Update the return value
    estim$lambda = resu_kcv$lambda
    estim$lambda.estimated = resu_kcv$lambda.estimated
    estim$lambda_var = resu_kcv$lambda_var
    estim$lambda_var.estimated = resu_kcv$lambda_var.estimated

    estim$Sigma = resu_kcv$Sigma
    estim$dof = resu_kcv$dof
    estim$dof.estimated = resu_kcv$dof.estimated

    myPsi = resu_kcv$Psi
    myPsi = resu_sbayes$Psi
    sigbar = ifelse(K>=2, mean(diag(resu_kcv$Sigma), na.rm = TRUE), resu_kcv$Sigma)
    estim$varresult <- convPsi2varresult(Psi = myPsi, Y = datY, X = datX,
                                         lambda = resu_kcv$lambda, scale_lambda = sigbar,
                                         type = type, ybar = ybar,
                                         Q_values = resu_kcv$q,
                                         callstr = cl
    )

    ## overwrite df.residual ##
    # Q_values = resu_kcv$q
    # sing_val = svd(sqrt(Q_values)*datX)$d #singular values
    # sigbar = ifelse(K>=2, mean(diag(resu_kcv$Sigma), na.rm = TRUE), resu_kcv$Sigma)
    # mykapp = sum((sing_val^2)/(sing_val^2 + resu_kcv$lambda/(1-resu_kcv$lambda)*sigbar*(N-1)/N),
    #              na.rm = TRUE) + as.vector(identical(type,'mean'))
    # for (i in 1:K) {
    #   estim$varresult[[i]]$df.residual = max(1, N - mykapp)
    # }

  }
  ##---------------------------------------------------

  ## Check return value ##
  ## components for 'varest'
  if (!with(estim, exists('varresult'))) {
    warning('VAR parameters were not estimated. Check the method.')
  }
  estim$datamat  = cbind(datY, datX)
  estim$y        = y
  estim$type     = type
  estim$p        = p
  estim$K        = K
  estim$obs      = N
  estim$totobs   = totobs
  estim$restrictions = NULL
  estim$call     = cl

  ## components for 'varshrinkest':
  estim$method   = method
  #$lambda
  #$lambda.estimated
  #$lambda_var
  #$lambda_var.estimated
  #$Sigma
  #$dof
  #$dof.estimated
  #...

  class(estim) <- c('varshrinkest', 'varest')

  return(estim)
}
