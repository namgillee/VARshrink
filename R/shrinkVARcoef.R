shrinkVARcoef <- function (Y, X, lambda, dof=Inf, prior_type='NCJ', TolDRes=1e-4, m0 = ncol(Y))
{
  # Shrinkage Estimation of VAR coefficients based on Bayesian estimation.
  #Likelihood model = multivariate t-distribution likelihood
  #Prior for coefficients = either conjugate normal distribution or non-conjugate normal distribution
  #prior for noise covariance = inverse Wishart distribution(s)
  #
  # Reference:
  #   Lee, N., Choi, H., Kim, S.-H., 2016.
  #   Bayes shrinkage estimation for high-dimensional VAR
  #   models with scale mixture of normal distributions for noise.
  #   Comput. Statist. Data Anal. 101, 250--276.
  #
  # 20 Nov. 2017, Namgil Lee, KNU

	#INPUTS
	#   X, Y : matrices
	#   lambda : a value between 0~1
  #   dof    : a value >= 1; If NULL or Inf, then use normal distribution
  #   prior_type : 'NCJ' for non-conjugate priors, 'CJ' for conjugate priors. Case insensitive.
  #   TolDRes
  #   m0     : a hyperparameter


  #Argument check
  if ((toupper(prior_type)!='NCJ')&&(toupper(prior_type)!='CJ')) {
    error('In:shrinkVARcoef:ArgChk','Incorrect argument:prior_type')
  }

  #Initialize
	MaxCount = 200 #iterative update of weight
	d = ncol(Y)
	K = ncol(X)
	p = K %/% d   #lag order
	N = nrow(Y)


  isNormal  <-  (is.null(dof)||is.infinite(dof))


  ## 1. When lambda==1, return zero matrix
	if (lambda==1) {
		Psihat = matrix(0, K, d)
		attr(Psihat, 'weight.estimated') <- FALSE
		attr(Psihat, 'weight') <- NULL
    attr(Psihat, 'noiseCov') <- NULL

    return(Psihat)
	}


  ## Initialize Psihat, Vhat, w via (weighted) least squares using Conjugate Priors
  w <- rep(1,N)  #Initial Weights
  l0 <- (m0+d+1)*rep(1,d) #InvWishart, L0 = diag(l0)
  if ((lambda <= 1e-7)&&(N<=(K+1))) {
    TMPs = svd( (t(X)%*%(X*w))/(N-1)*(1-lambda) + lambda*diag(rep(1,K)) )
    TMPY = (t(X)%*%(Y*w))/(N-1)*(1-lambda)
    Psihat = TMPs$v %*% ((t(TMPs$u) %*% TMPY)/TMPs$d)
  } else {
    U <- chol( (t(X)%*%(X*w))/(N-1)*(1-lambda) + lambda*diag(rep(1,K)) )
    Psihat <- backsolve(U, backsolve(U, (t(X)%*%(Y*w))/(N-1)*(1-lambda), transpose=TRUE) )
  }
  Vhat <- t(Y) %*% ((Y - (X%*%Psihat))*w)
  diag(Vhat) <- diag(Vhat) + l0
  Vhat <- Vhat / (m0+N+d+1)


  ## 2. Normal distribution & Conjugate Prior: No iteration
  if (isNormal && (toupper(prior_type)=='CJ'))  {
		attr(Psihat, 'weight.estimated') <- FALSE
		attr(Psihat, 'weight') <- w
		attr(Psihat, 'noiseCov') <- Vhat
    return(Psihat)
	}


  ## 3. Iteration for w, Psihat, Vhat
  if (isNormal) {
    h <- function(x) 1
    flag_for_ncj_loop <- TRUE #If isNormal && 'NCJ', run ncj without update of weight
  } else {
    h <- function(x) (dof + d)/(dof + x)
    flag_for_ncj_loop <- FALSE #If !isNormal, run cj (initially) with update of weight
  }
  eVhat <- eigen(Vhat)


  #AT LEAST ONE MORE UPDATE OF w,Psihat,Vhat IS NEEDED
  for (count in 1:MaxCount) {
    iVhat <- eVhat$vectors %*% diag(1/pmax((eVhat$values),1e-18)) %*% t(eVhat$vectors)

		### w is determined by w_t = h( e_t' V^{-1} e_t )
    if (!isNormal && !flag_for_ncj_loop) {
      w_prev     <- w
      Err        <- Y - X%*%Psihat
      w          <- h( colSums(t(Err) * (iVhat %*% t(Err))) )
    }

    ### Psihat
    if (toupper(prior_type)=='CJ' || !flag_for_ncj_loop) {
      if ((lambda <= 1e-7)&&(N<=(K+1))) {
        TMPs = svd( (t(X)%*%(X*w))/(N-1)*(1-lambda) + lambda*diag(rep(1,K)) )
        TMPY = (t(X)%*%(Y*w))/(N-1)*(1-lambda)
        Psihat = TMPs$v %*% ((t(TMPs$u) %*% TMPY)/TMPs$d)
      } else {
        U <- chol( (t(X)%*%(X*w))/(N-1)*(1-lambda) + lambda*diag(rep(1,K)) )
        Psihat <- backsolve(U, backsolve(U, (t(X)%*%(Y*w))/(N-1)*(1-lambda), transpose=TRUE) )
      }

    } else { #'NCJ' && flag_for_ncj_loop
      if ( !exists('eXX') || (!isNormal && !flag_for_ncj_loop)) {
        eXX <- eigen(t(X)%*%(X*w))
      }
  		mXY <- t(X)%*%(Y*w)%*%iVhat
  		Dv = 1/( kronecker(1/pmax((eVhat$values),1e-18) , eXX$values) + lambda/(1-lambda)*(N-1) )
  		Psihat = matrix(Dv,K,d) * ( t(eXX$vectors) %*% mXY %*% eVhat$vectors  )
  		Psihat = eXX$vectors %*% (Psihat %*% t(eVhat$vectors))
		}

    ### Vhat
		ev_prev <- eVhat$values
		Vhat <- t(Y) %*% ((Y - (X%*%Psihat))*w)
		diag(Vhat) <- diag(Vhat) + l0
		Vhat <- Vhat / (m0+N+d+1)
		Vhat <- (Vhat+t(Vhat))/2   #To ensure real eigenvalues
		eVhat <- eigen(Vhat)


    #Convergence Criterion
    if (isNormal || flag_for_ncj_loop) { #w is not changed, so we use eVhat$values
      if ( sum(abs((ev_prev - eVhat$values)^2)) <= TolDRes * sum(ev_prev^2) ) {
        break
      }

    } else {
      if ( sum((w - w_prev)^2) <= TolDRes^2 * sum(w_prev^2) ) {
          if (toupper(prior_type)=='CJ') {
            break
          } else {
            flag_for_ncj_loop = TRUE
          }
      }
    }
  }
  attr(Psihat, 'weight.estimated') <- TRUE
  attr(Psihat, 'weight') <- w
  attr(Psihat, 'noiseCov') <- Vhat

  return(Psihat)

}
