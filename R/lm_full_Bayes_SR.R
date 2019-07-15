lm_full_Bayes_SR <- function(Y, X, dof=Inf, burnincycle=1000, mcmccycle=2000)
{
  #lm_full_Bayes_SR <- function(Y, X, dof=Inf, burnincycle=1000, mcmccycle=2000)
  # Gibbs MCMC algorithm for sampling and estimating VAR parameters from 
  # conditional posterior distributions. 
  #
  # The likelihood function (distribution for noise vector) 
  #is multivariate t-distribution with degree of freedom dof=nu \geq 0.
  #
  # The priors are noninformative priors: 
  #   1) the shrinkage prior for the VAR coefficients, and 
  #   2) the reference prior for the VAR noise covariance matrix. 
  #
  # For further detail, see Ni & Sun (2005). 
  #
  # Note: Lag order p cannot be estimated but fixed, and 
  #       it is inferred from data matrix Y and X (see below)
  #
  # INPUTS: 
  #   Y, X: Y is an N x D matrix, 
  #         X is an N x K matrix, 
  #         The size K is either (K = D^2p, or, K = D(Dp+1)). 
  #         So, p can be estimated by 
  #            p == K %/% D
  #
  #   dof : If dof = Inf, then we apply multivarate normal distribution, and do not estimate Q. 
  #         If dof = 'a finite positive value', then we apply the multivariate t-distribution 
  #                  with fixed dof, and estimate Q
  #         If dof = NULL or dof <= 0, we estimate both dof and Q
  #
  # OUTPUT
  #   Returns a list of two VARparam objects, 
  #   one for estimates of Psi, Sig, dof, etc., 
  #   and the other one for SE of the estimates.
  #   See Ni & Sun (2005) for the posterior mean estimate 
  #   and LINEXloss estimate.
  #   
  #   $param$Psi :   The coefficient matrix Psi computed by a posterior mean
  #   $param$Sigma : The Sigma matrix computed as a posteror mean
  #   $param$dof : The dof value, which is either estimated or fixed
  #   $param$delta
  #   $param$lambda
  #   $se.param$Psi, 
  #   $se.param$Sigma, 
  #   $se.param$dof, 
  #   $se.param$delta, 
  #   $se.param$lambda
  #           : The standard errors for the posterior mean estimates, 
  #             which is computed by sd_of_samples / sqrt(mcmccycle). 
  #             where sd_of_samples is the sd of MCMC samples. 
  #             If dof is fixed (given as an argument), then 
  #             se.param$dof==0. 
  #   $LINEXVARmodel : list with $Psi and $SE.Psi, 
  #           representing a LINEX loss estimate with aij=-4. 
  #           See, loss #2 in Ni & Sun (2005)
  #
  # Last modified at Nov. 8th, 2018 by Namgil Lee, Kangwon National University, South Korea.
  # Reference: Ni and Sun (2005)
  
  # Find the lag order p
  d = ncol(Y)   #size D
  K = ncol(X)   #size K = D^2p, D(Dp+1)
  p = K %/% d   #lag order
  
  N = nrow(X)
	J = d*K

	#initial values
  delta = 0.001
  
  estimate_dof = FALSE
  if (is.null(dof)) {
    estimate_dof = TRUE
    require(ars)
    w = Inf
  } else if (dof<=0) {
    estimate_dof = TRUE
    dof = NULL
    w = Inf
  } else {
	  w = dof/2 #5 #w=nu/2
  }
	Q = rep(1,N) #weight vector
	vec_a4 <- rep(-4,J)
	vec_a4[1+(1+p*d)*(0:(d-1))] <- 0.001
	#vec_a8 <- rep(-8,J)##################not used##
	#vec_a8[1+(1+p*d)*(0:(d-1))] <- 0.001#not used##
	Sig = 0.0001 + diag(1/rgamma(d,shape=3,scale=1/4), d) # inverse gamma

  #Gibbs sampler: burn-in cycles & MCMC cycles
	Psi01 <- Psi02 <- Sig01 <- dof01 <- delta01 <- 0
	Psi01SE <- Psi02SE <- Sig01SE <- dof01SE <- delta01SE <- 0
	
	for (i in 1:(burnincycle+mcmccycle))
	{
	  #update phi, delta and Sig:

    ##### (1) draw phi from a multivariate normal distribution #####
	  # mu_Q = delta * (Sig \otimes (X'QX)^{-1} + delta I_J)^{-1} vec(hat{phi}_Q)
	  #            = ( (1/delta) I \otimes I  +  Sig^{-1} \otimes X'QX )^{-1} vec(X'QY Sig^{-1})
	  #            = V_Q vec(X'QY Sig^{-1})
	  # phi ~ N_J(mu_Q, V_Q)
	  # Suppose that V_Q = V = UDU', then, U = kron(Us, Ux)
	  #	U'phi ~ N_J(U' mu, D)
	  #               = N_J(U'V vecXY, D)
	  #               = N_J(DU' vecXY, D)
	  #    Since U = kron(Us, Ux), we have 
	  #    U'phi ~ N_J (D vec(Ux' (X'QY Sig^{-1}) Us), D)
	  #
	  if (i==1) {
	  	eSig <- eigen(Sig)						# eigenvalue decomposition of Sig
		  invSig <- eSig$vectors %*% diag(1/pmax(eSig$values,1e-20)) %*% t(eSig$vectors)
	  }
	  eXX <- eigen(t(X)%*%(X*Q))			# eigenvalue decomposition of X'QX
	  mXY <- t(X)%*%(Y*Q)%*%invSig		# matrix X'QY Sig^{-1}

	  Dv = 1/( kronecker(1/pmax(eSig$values,1e-14) , pmax(eXX$values,1e-14)) + 1/delta) # eigenvalues D=Dv of V_Q
	  #Dv = 1/( kronecker(1/abs(eSig$values) , abs(eXX$values)) + 1/delta) # eigenvalues D=Dv of V_Q
	  
	  mu0 = Dv * as.vector( t(eXX$vectors) %*% mXY %*% eSig$vectors  ) # mean in U'phi ~ N_J(mu0, Dv)
	  phi = mu0 + sqrt(Dv) * rnorm(J, 0, 1)		# U'phi
	  dim(phi) = c(d*p+1, d)
	  phi = eXX$vectors %*% (phi %*% t(eSig$vectors))  # U * U'phi
	  ######################################


	  ##### (2) draw delta from InvGamma(J/2-1, x/2) #####
	  x <- sum(phi^2)
	  if (x < 1e-20)
	    delta = max(x/J, 1e-30)
	  else
      delta = 1/rgamma(1,shape=J/2-1,scale=2/x) # inverse gamma
	  ######################################


    ##### (3) draw Sig based on Sun & Ni (2004) #####
    Sk <- Y - X%*%phi
	  Sk <- t(Sk)%*%(Sk*Q)					# S_k

	  logLambda <- log(pmax(eSig$values,1e-20)) #eSig <- eigen(Sig)
	  SigStar <- (eSig$vectors)%*%diag(logLambda)%*%t(eSig$vectors) # Sig_star

          zij <- rnorm(d*(d+1)/2)					# z_ij are standard normals.
          zij <- zij/sqrt(sum(zij^2))				# v_ij are upper-triangular part of V, 
	  V <- matrix(0, d, d)						#that are normalized standard normals.
	  V[upper.tri(V,diag=TRUE)] <- zij			# 
	  for (j in 2:d)							#
	    for (k in 1:(j-1))						#
	      V[j,k] <- V[k,j]						# V is a symmetric metric

    W <- SigStar + rnorm(1)*V				# W = Sigstar + z*V
    eW <- eigen(W)						
	  id_sort <- sort.list(eW$values,decreasing=TRUE)
	  eW$values <- eW$values[id_sort]
	  eW$vectors <- eW$vectors[,id_sort]
	  Cstar <- eW$values 

	  alpha_k <- N/2*sum(logLambda-Cstar) + 
	      1/2*sum( (invSig - eW$vectors%*%diag(1/exp(Cstar))%*%t(eW$vectors)) * Sk )
	  for (j in 1:(d-1)) {
	    for (k in (j+1):d) {
                alpha_k <- alpha_k + log(logLambda[j]-logLambda[k]) - log(Cstar[j]-Cstar[k])
		if (is.nan(alpha_k))
		  alpha_k = -Inf
	    }
	  }

	  if ( runif(1) <= min(1,exp(alpha_k)) ) {		# Update Sig, eSig, invSig
      Sig <- eW$vectors%*%diag(exp(Cstar))%*%t(eW$vectors)
    	invSig <- eW$vectors%*%diag( 1/pmax(exp(Cstar),1e-20) )%*%t(eW$vectors)
	  	eSig <- list(vectors=eW$vectors, values=exp(Cstar))
	  } #if not, use previous eSig 
  	######################################

	  ##### (4) draw Q from gamma #####
	  if (is.null(dof) || !is.infinite(dof) )
	  {
	  	#If dof is Inf, then it is a multivarate normal distribution, and do not update Q. 
	  	#If dof is not Inf, then it is a multivariate t-distribution, and update Q
	    if (is.infinite(w)) {
	      Q = rep(1,N) 
	    } else {
	      x = (Y - X %*% phi)
	      x = w + 0.5* rowSums(x * (x %*% invSig))				  # x : rate, beta, 1/scale
	      Q = rgamma(N, shape=(w + 0.5*d), scale=1) / x	#alpha=(0.5(nu + d)), scale = 1/x  
	    }
	  }
		######################################

    ##### (5) draw w by MCMC #####
    if (estimate_dof) { ##if(is.null(dof))
      f_log = function(x,N,Q) {
	  N*x*log(x) + x*sum(log(Q),na.rm=TRUE) - N*lgamma(x) - (1+sum(Q))*x
        #N*x*log(x) + x*log(prod(Q)) - N*lgamma(x) - (1+sum(Q))*x
        #N*exp(x)*x + exp(x)*log(prod(Q)) - N*lgamma(exp(x)) - (1+sum(Q))*exp(x) + x
      }
      f_dif = function(x,N,Q) {
        N*log(x) + N + sum(log(Q),na.rm=TRUE) - N*digamma(x) - (1+sum(Q))
        #N*log(x) + N + log(prod(Q)) - N*digamma(x) - (1+sum(Q))
        #(N*x + N + log(prod(Q)) - N*digamma(exp(x)) - (1+sum(Q))) * exp(x) + 1
      }
      tmp = capture.output( {
             w = ars(n=1, f=f_log, fprima=f_dif, x=10, m=1, lb=TRUE, xlb=1e-15, N=N, Q=Q)
      })
      #w = ars(n=1, f=f_log, fprima=f_dif, x=2, m=1, N=N, Q=Q)
      #w = exp(w)
      #ub=TRUE, xub=1e15
    }
    ##############################

	  ##### update Psi01, Psi02, Sig01 #####
    if (i>=(burnincycle+1)) {
      
        Psi01 <- Psi01 + phi / mcmccycle
	      Psi02 <- Psi02 + exp(-vec_a4 * phi) / mcmccycle
	      Sig01 <- Sig01 + Sig / mcmccycle
        delta01 <- delta01 + delta / mcmccycle
        dof01 <- dof01 + w*2/mcmccycle
        
        # Note that, SE^2 = S^2/n = sum(x_i^2/n/(n-1)) - xbar^2/(n-1)
        Psi01SE <- Psi01SE + phi^2 / mcmccycle / (mcmccycle-1)
        Psi02SE <- Psi02SE + (exp(-vec_a4 * phi))^2 / mcmccycle / (mcmccycle-1)
        Sig01SE <- Sig01SE + Sig^2 / mcmccycle / (mcmccycle-1)
        delta01SE <- delta01SE + delta^2 / mcmccycle / (mcmccycle-1)
        if (estimate_dof)
          dof01SE <- dof01SE + (w*2)^2 /mcmccycle / (mcmccycle-1)
        
	  }
	  ##########

	}

	#Compute Psi02 first
	mySE.Psi02 = matrix( sqrt( Psi02SE - Psi02^2 / (mcmccycle-1) ) ,K,d) #SE value (before -log)
	Psi02 = -log(Psi02)/vec_a4     #(after -log)
	myPsi02 = matrix(Psi02, K, d)  #(after -log)
	mySE.Psi02 = mySE.Psi02 / abs(vec_a4)  #Lipshitz constant: 1/a
	
	
	#######
	
	# Collect return values
	myPsi01 = matrix(Psi01, K, d)   # Psi01 matrix
	mySE.Psi01 = matrix( sqrt( Psi01SE - Psi01^2 / (mcmccycle-1) ) ,K,d)
	mySigma = Sig01
  mySE.Sigma = sqrt( Sig01SE - Sig01^2 / (mcmccycle-1) )
  if (estimate_dof) {
    mydof = dof01
    mySE.dof = sqrt( dof01SE - dof01^2 / (mcmccycle-1) )
  } else {
    mydof = dof
    mySE.dof = 0
  }
  # Estimate Q values by the mode of posterior, (shape-1)*scale (see, lines 175--186)
  if (is.infinite(mydof)) {
    myq = rep(1,N)
  } else {
    x = (Y - X %*% myPsi01)
    x = 0.5*mydof + 0.5*rowSums(x * t(solve( mySigma, t(x) )))
    myq = (0.5*mydof + 0.5*d - 1)/min(x, 1e10) 
  }

  res <- NULL
  res$Psi = myPsi01
  res$Sigma = mySigma
  res$dof = mydof
  res$delta = delta01
  res$lambda = 1/delta01
  res$lambda.estimated = TRUE
  res$dof.estimated = estimate_dof
  res$q = myq
  
  res$se.param <- NULL
  res$se.param$Psi = mySE.Psi01
  res$se.param$Sigma = mySE.Sigma
  res$se.param$dof = mySE.dof
  res$se.param$delta = sqrt( delta01SE - delta01^2 / (mcmccycle-1) )
  res$se.param$lambda = res$se.param$delta / delta01^2
 
	res$LINEXparam <- NULL
	res$LINEXparam$Psi = myPsi02
	res$se.LINEXparam <- NULL
	res$se.LINEXparam$Psi = mySE.Psi02
	
	res
}