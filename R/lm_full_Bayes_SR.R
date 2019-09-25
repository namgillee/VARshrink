#' Full Bayesian Shrinkage Estimation Method for Multivariate Regression
#'
#' Estimate regression coefficients and scale matrix for noise by using
#' Gibbs MCMC algorithm. The function assumes 1) multivariate t-distribution
#' for noise as a sampling distribution, and 2) noninformative priors for
#' regression coefficients and scale matrix for noise.
#'
#' Consider the multivariate regression:
#' \deqn{Y = X Psi + e, \quad e ~ mvt(0, dof, Sigma).}
#' Psi is a M-by-K matrix of regression coefficients and
#' Sigma is a K-by-K scale matrix for multivariate t-distribution for noise.
#'
#' Sampling distribution for noise e is multivariate t-distribution with
#' degree of freedom dof and scale matrix Sigma: e ~ mvt(0, dof, Sigma).
#' The priors are noninformative priors: 1) the shrinkage prior for regression
#' coefficients Psi, and 2) the reference prior for scale matrix Sigma.
#'
#' The function implements Gibbs MCMC algorithm for estimating regression
#' coefficients Psi and scale matrix Sigma.
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param dof Degree of freedom for multivariate t-distribution.
#' If dof = Inf (default), then multivariate normal distribution is applied and
#' weight vector q is not estimated. If dof = NULL or dof <= 0, then dof and q
#' are estimated automatically. If dof is a positive number, q is estimated.
#' @param burnincycle,mcmccycle Number of burnin cycles is the number of
#' initially generated sample values to drop. Number of MCMC cycles is the
#' number of generated sample values to compute estimates.
#' @return A list object with estimated parameters: Psi, Sigma, dof, delta
#' (delta is the reciprocal of lambda), and lambda.
#' Additional components are se.param (standard error of the parameters) and
#' LINEXVARmodel (estimates under LINEX loss).
#' @references S. Ni and D. Sun (2005). Bayesian estimates for vector
#' autoregressive models. Journal of Business & Economic Statistics 23(1),
#' 105-117.
#' @importFrom stats rgamma rnorm runif
#' @importFrom utils capture.output
#
# Last modified: Nov. 8th, 2018 by Namgil Lee @ Kangwon National University,
# South Korea.

lm_full_Bayes_SR <- function(Y, X, dof = Inf, burnincycle = 1000,
                             mcmccycle = 2000) {
  # Find the lag order p
  K <- ncol(Y)   #size D
  M <- ncol(X)   #size M = D^2p, D(Dp+1)
  p <- M %/% K   #lag order

  N <- nrow(X)
	J <- K * M

	#initial values
  delta <- 0.001

  estimate_dof <- FALSE
  if (is.null(dof)) {
    estimate_dof <- TRUE
    w <- Inf
  } else if (dof <= 0) {
    estimate_dof <- TRUE
    dof <- NULL
    w <- Inf
  } else {
	  w <- dof / 2
  }
	Q <- rep(1, N) #weight vector
	vec_a4 <- rep(-4, J)
	vec_a4[1 + M * (0:(K - 1))] <- 0.001
	Sig <- 0.0001 + diag(1 / rgamma(K, shape = 3, scale = 1 / 4), K) #inverse gamma

  #Gibbs sampler: burn-in cycles & MCMC cycles
	Psi01 <- Psi02 <- Sig01 <- dof01 <- delta01 <- 0
	Psi01SE <- Psi02SE <- Sig01SE <- dof01SE <- delta01SE <- 0

	for (i in 1:(burnincycle + mcmccycle)) {
	  #update phi, delta and Sig:

    ##### (1) draw phi from a multivariate normal distribution #####
	  # mu_Q = delta * (Sig \otimes (X'QX)^{-1} + delta I_J)^{-1} vec(hat{phi}_Q)
	  #            = ( (1/delta) I \otimes I  +  Sig^{-1} \otimes X'QX )^{-1} %*%
	  #              vec(X'QY Sig^{-1})
	  #            = V_Q vec(X'QY Sig^{-1})
	  # phi ~ N_J(mu_Q, V_Q)
	  # Suppose that V_Q = V = UDU', then, U = kron(Us, Ux)
	  #	U'phi ~ N_J(U' mu, D)
	  #               = N_J(U'V vecXY, D)
	  #               = N_J(DU' vecXY, D)
	  #    Since U = kron(Us, Ux), we have
	  #    U'phi ~ N_J (D vec(Ux' (X'QY Sig^{-1}) Us), D)
	  #
	  if (i == 1) {
	  	eSig <- eigen(Sig)						# eigenvalue decomposition of Sig
		  invSig <- eSig$vectors %*% diag(1 / pmax(eSig$values, 1e-20)) %*%
		    t(eSig$vectors)
	  }
	  eXX <- eigen(t(X) %*% (X * Q))			# eigenvalue decomposition of X'QX
	  mXY <- t(X) %*% (Y * Q) %*% invSig	# matrix X'QY Sig^{-1}

	  Dv <- 1 / (kronecker(1 / pmax(eSig$values, 1e-14),
	                       pmax(eXX$values, 1e-14)) + 1 / delta)
	  # eigenvalues D=Dv of V_Q
	  #Dv = 1/( kronecker(1/abs(eSig$values) , abs(eXX$values)) + 1/delta)
	  # eigenvalues D=Dv of V_Q

	  mu0 <- Dv * as.vector(t(eXX$vectors) %*% mXY %*% eSig$vectors)
	  # mean in U'phi ~ N_J(mu0, Dv)
	  phi <- mu0 + sqrt(Dv) * rnorm(J, 0, 1)		# U'phi
	  dim(phi) <- c(K * p + 1, K)
	  phi <- eXX$vectors %*% (phi %*% t(eSig$vectors))  # U * U'phi
	  ######################################


	  ##### (2) draw delta from InvGamma(J/2-1, x/2) #####
	  x <- sum(phi^2)
	  if (x < 1e-20)
	    delta <- max(x / J, 1e-30)
	  else
      delta <- 1 / rgamma(1, shape = J / 2 - 1, scale = 2 / x) # inverse gamma
	  ######################################


    ##### (3) draw Sig based on Sun & Ni (2004) #####
    Sk <- Y - X %*% phi
	  Sk <- t(Sk) %*% (Sk*Q)					# S_k

	  logLambda <- log(pmax(eSig$values, 1e-20)) #eSig <- eigen(Sig)
	  SigStar <- (eSig$vectors) %*% diag(logLambda) %*% t(eSig$vectors) # Sig_star

    zij <- rnorm(K * (K + 1) / 2)					# z_ij are standard normals.
    zij <- zij / sqrt(sum(zij^2))				# v_ij are upper-triangular part of V,
	  V <- matrix(0, K, K)						#that are normalized standard normals.
	  V[upper.tri(V, diag = TRUE)] <- zij			#
	  for (j in 2:K)							#
	    for (k in 1:(j - 1))						#
	      V[j, k] <- V[k, j]						# V is a symmetric metric

    W <- SigStar + rnorm(1) * V				# W = Sigstar + z*V
    eW <- eigen(W)
	  id_sort <- sort.list(eW$values, decreasing = TRUE)
	  eW$values <- eW$values[id_sort]
	  eW$vectors <- eW$vectors[, id_sort]
	  Cstar <- eW$values

	  alpha_k <- N / 2 * sum(logLambda - Cstar) +
      1/2*sum((invSig - eW$vectors %*% diag(1 / exp(Cstar)) %*% t(eW$vectors)) * Sk)
	  for (j in 1:(K - 1)) {
	    for (k in (j + 1):K) {
        alpha_k <- alpha_k + log(logLambda[j] - logLambda[k]) -
          log(Cstar[j] - Cstar[k])
        if (is.nan(alpha_k))
          alpha_k <- -Inf
	    }
	  }

	  if (runif(1) <= min(1, exp(alpha_k))) {
	    # Update Sig, eSig, invSig
      Sig <- eW$vectors %*% diag(exp(Cstar)) %*% t(eW$vectors)
    	invSig <- eW$vectors %*% diag(1 / pmax(exp(Cstar), 1e-20)) %*%
    	  t(eW$vectors)
	  	eSig <- list(vectors = eW$vectors, values = exp(Cstar))
	  } #if not, use previous eSig
  	######################################

	  ##### (4) draw Q from gamma #####
	  if (is.null(dof) || !is.infinite(dof))
	  {
	  	#If dof is Inf, then it is a multivarate normal distribution,
	    #and do not update Q.
	  	#If dof is not Inf, then it is a multivariate t-distribution, and update Q
	    if (is.infinite(w)) {
	      Q <- rep(1, N)
	    } else {
	      x <- (Y - X %*% phi)
	      x <- w + 0.5 * rowSums(x * (x %*% invSig))			 # x : rate, beta, 1/scale

	      #gamma distribution with alpha=(0.5(nu + K)), scale = 1/x
	      Q <- rgamma(N, shape = (w + 0.5 * K), scale = 1) / x
	    }
	  }
		######################################

    ##### (5) draw w by MCMC #####
    if (estimate_dof) {
      f_log <- function(x, N, Q) {
        N * x * log(x) + x * sum(log(Q), na.rm = TRUE) - N * lgamma(x) -
          (1 + sum(Q)) * x
        #N*x*log(x) + x*log(prod(Q)) - N*lgamma(x) - (1+sum(Q))*x
        #N*exp(x)*x + exp(x)*log(prod(Q)) - N*lgamma(exp(x)) -
        #  (1+sum(Q))*exp(x) + x
      }
      f_dif <- function(x, N, Q) {
        N * log(x) + N + sum(log(Q), na.rm = TRUE) - N * digamma(x) -
          (1 + sum(Q))
        #N*log(x) + N + log(prod(Q)) - N*digamma(x) - (1+sum(Q))
        #(N*x + N + log(prod(Q)) - N*digamma(exp(x)) - (1+sum(Q))) * exp(x) + 1
      }
      tmp <- capture.output( {
        w <- ars::ars(n = 1, f = f_log, fprima = f_dif, x = 10, m = 1,
                lb = TRUE, xlb = 1e-15, N = N, Q = Q)
      })
      #w = ars::ars(n=1, f=f_log, fprima=f_dif, x=2, m=1, N=N, Q=Q)
      #w = exp(w)
      #ub=TRUE, xub=1e15
    }
    ##############################

	  ##### update Psi01, Psi02, Sig01 #####
    if (i >= (burnincycle + 1)) {

        Psi01 <- Psi01 + phi / mcmccycle
	      Psi02 <- Psi02 + exp(-vec_a4 * phi) / mcmccycle
	      Sig01 <- Sig01 + Sig / mcmccycle
        delta01 <- delta01 + delta / mcmccycle
        dof01 <- dof01 + w * 2 / mcmccycle

        # Note that, SE^2 = S^2/n = sum(x_i^2/n/(n-1)) - xbar^2/(n-1)
        Psi01SE <- Psi01SE + phi^2 / mcmccycle / (mcmccycle - 1)
        Psi02SE <- Psi02SE + (exp(-vec_a4 * phi))^2 / mcmccycle / (mcmccycle - 1)
        Sig01SE <- Sig01SE + Sig^2 / mcmccycle / (mcmccycle - 1)
        delta01SE <- delta01SE + delta^2 / mcmccycle / (mcmccycle - 1)
        if (estimate_dof)
          dof01SE <- dof01SE + (w * 2)^2 / mcmccycle / (mcmccycle - 1)

	  }
	  ##########

	}

	#Compute Psi02 first
	#mySE.Psi02: SE value (before -log)
	mySE.Psi02 <- matrix(sqrt(Psi02SE - Psi02^2 / (mcmccycle - 1)), M, K)
	Psi02 <- -log(Psi02)/vec_a4     #(after -log)
	myPsi02 <- matrix(Psi02, M, K)  #(after -log)
	mySE.Psi02 <- mySE.Psi02 / abs(vec_a4)  #Lipshitz constant: 1/a

	#######

	# Collect return values
	myPsi01 <- matrix(Psi01, M, K)   # Psi01 matrix
	mySE.Psi01 <- matrix(sqrt(Psi01SE - Psi01^2 / (mcmccycle - 1)), M, K)
	mySigma <- Sig01
  mySE.Sigma <- sqrt(Sig01SE - Sig01^2 / (mcmccycle - 1))
  if (estimate_dof) {
    mydof <- dof01
    mySE.dof <- sqrt(dof01SE - dof01^2 / (mcmccycle - 1))
  } else {
    mydof <- dof
    mySE.dof <- 0
  }
  # Estimate Q values by the mode of posterior, (shape-1)*scale
  #(see, lines 175--186)
  if (is.infinite(mydof)) {
    myq <- rep(1, N)
  } else {
    x <- (Y - X %*% myPsi01)
    x <- 0.5 * mydof + 0.5 * rowSums(x * t(solve(mySigma, t(x))))
    myq <- (0.5 * mydof + 0.5 * K - 1) / min(x, 1e10)
  }

  res <- NULL
  res$Psi <- myPsi01
  res$Sigma <- mySigma
  res$dof <- mydof
  res$delta <- delta01
  res$lambda <- 1 / delta01
  res$lambda.estimated <- TRUE
  res$dof.estimated <- estimate_dof
  res$q <- myq

  res$se.param <- NULL
  res$se.param$Psi <- mySE.Psi01
  res$se.param$Sigma <- mySE.Sigma
  res$se.param$dof <- mySE.dof
  res$se.param$delta <- sqrt(delta01SE - delta01^2 / (mcmccycle - 1))
  res$se.param$lambda <- res$se.param$delta / delta01^2

	res$LINEXparam <- NULL
	res$LINEXparam$Psi <- myPsi02
	res$se.LINEXparam <- NULL
	res$se.LINEXparam$Psi <- mySE.Psi02

	res
}
