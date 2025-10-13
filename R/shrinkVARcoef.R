#' Semiparametric Bayesian Shrinkage Estimator for
#' Multivariate Regression
#'
#' Compute the semiparametric Bayesian shrinkage estimator of Psi and Sigma
#' for a given shrinkage parameter lambda.
#' The function is a private function for lm_semi_Bayes_PCV() and
#' lm_ShVAR_KCV().
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param lambda A shrinkage intensity parameter value between 0~1.
#' @param dof Degree of freedom for multivariate t-distribution.
#' If NULL or Inf, then use multivariate normal distribution.
#' @param prior_type "NCJ" for non-conjugate prior and "CJ" for conjugate
#' prior for scale matrix Sigma.
#' @param TolDRes Tolerance parameter for stopping criterion.
#' @param m0 A hyperparameter for inverse Wishart distribution for Sigma
#' @references N. Lee, H. Choi, and S.-H. Kim (2016). Bayes shrinkage
#' estimation for high-dimensional VAR models with scale mixture of normal
#' distributions for noise. Computational Statistics & Data Analysis 101,
#' 250-276. doi: 10.1016/j.csda.2016.03.007
#
# Last modified: 20 Nov. 2017, Namgil Lee @ Kangwon National University

shrinkVARcoef <- function(Y, X, lambda, dof = Inf, prior_type = "NCJ",
                          TolDRes = 1e-4, m0 = ncol(Y)) {
  #Argument check
  if ((toupper(prior_type) != "NCJ") && (toupper(prior_type) != "CJ")) {
    stop(paste("Unknown argument prior_type =", prior_type))
  }

  #Initialize
	MaxCount <- 200 #iterative update of weight
	d <- ncol(Y)
	K <- ncol(X)
	N <- nrow(Y)

  isNormal  <-  (is.null(dof) || is.infinite(dof))

  ## 1. When lambda==1, return zero matrix
	if (lambda == 1) {
    Psihat <- matrix(0, K, d)
		attr(Psihat, "weight.estimated") <- FALSE
		attr(Psihat, "weight") <- NULL
    attr(Psihat, "noiseCov") <- NULL

    return(Psihat)
	}


  ## Initialize Psihat, Vhat, w via (weighted) least squares using
  ## Conjugate Priors
  w <- rep(1, N)  #Initial Weights
  l0 <- (m0 + d + 1) * rep(1, d) #InvWishart, L0 = diag(l0)
  if ((lambda <= 1e-7) && (N <= (K + 1))) {
    TMPs <- svd((t(X) %*% (X * w)) / (N - 1) * (1 - lambda) +
                  lambda * diag(rep(1, K)))
    TMPY <- (t(X) %*% (Y * w)) / (N - 1) * (1 - lambda)
    Psihat <- TMPs$v %*% ((t(TMPs$u) %*% TMPY) / TMPs$d)
  } else {
    U <- chol((t(X) %*% (X * w)) / (N - 1) * (1 - lambda) +
                lambda * diag(rep(1, K)))
    Psihat <- 
      backsolve(U,
                backsolve(U, (t(X) %*% (Y * w)) / (N - 1) *
                            (1 - lambda), transpose = TRUE))
  }
  Vhat <- t(Y) %*% ((Y - (X %*% Psihat)) * w)
  diag(Vhat) <- diag(Vhat) + l0
  Vhat <- Vhat / (m0 + N + d + 1)


  ## 2. Normal distribution & Conjugate Prior: No iteration
  if (isNormal && (toupper(prior_type) == "CJ"))  {
		attr(Psihat, "weight.estimated") <- FALSE
		attr(Psihat, "weight") <- w
		attr(Psihat, "noiseCov") <- Vhat
    return(Psihat)
	}


  ## 3. Iteration for w, Psihat, Vhat
  if (isNormal) {
    #If isNormal && 'NCJ', run ncj without update of weight
    h <- function(x) 1
    flag_for_ncj_loop <- TRUE
  } else {
    #If !isNormal, run cj (initially) with update of weight
    h <- function(x) (dof + d) / (dof + x)
    flag_for_ncj_loop <- FALSE
  }
  eVhat <- eigen(Vhat)


  #AT LEAST ONE MORE UPDATE OF w,Psihat,Vhat IS NEEDED
  for (count in 1:MaxCount) {
    iVhat <- eVhat$vectors %*% diag(1 / pmax((eVhat$values), 1e-18)) %*%
      t(eVhat$vectors)

    ### w is determined by w_t = h( e_t' V^{-1} e_t )
    if (!isNormal && !flag_for_ncj_loop) {
      w_prev     <- w
      Err        <- Y - X %*% Psihat
      w          <- h(colSums(t(Err) * (iVhat %*% t(Err))))
    }

    ### Psihat
    if (toupper(prior_type) == "CJ" || !flag_for_ncj_loop) {
      if ((lambda <= 1e-7) && (N <= (K + 1))) {
        TMPs <- svd((t(X) %*% (X * w)) / (N - 1) * (1 - lambda) +
                      lambda * diag(rep(1, K)))
        TMPY <- (t(X) %*% (Y * w)) / (N - 1) * (1 - lambda)
        Psihat <- TMPs$v %*% ((t(TMPs$u) %*% TMPY) / TMPs$d)
      } else {
        U <- chol((t(X) %*% (X * w)) / (N - 1) * (1 - lambda) +
                    lambda * diag(rep(1, K)))
        Psihat <- backsolve(U, backsolve(U, (t(X) %*% (Y * w)) / (N - 1) *
                                           (1 - lambda), transpose = TRUE))
      }

    } else { #'NCJ' and flag_for_ncj_loop
      if (!exists("eXX") || (!isNormal && !flag_for_ncj_loop)) {
        eXX <- eigen(t(X) %*% (X * w))
      }
      mXY <- t(X) %*% (Y * w) %*% iVhat
  		Dv <- 1 / (kronecker(1 / pmax((eVhat$values), 1e-18), eXX$values) +
  		           lambda / (1 - lambda) * (N - 1))
      Psihat <- matrix(Dv, K, d) * (t(eXX$vectors) %*% mXY %*% eVhat$vectors)
      Psihat <- eXX$vectors %*% (Psihat %*% t(eVhat$vectors))
		}

    ### Vhat
    ev_prev <- eVhat$values
		Vhat <- t(Y) %*% ((Y - (X %*% Psihat)) * w)
		diag(Vhat) <- diag(Vhat) + l0
		Vhat <- Vhat / (m0 + N + d + 1)
		Vhat <- (Vhat + t(Vhat)) / 2 #To ensure real eigenvalues
		eVhat <- eigen(Vhat)

    # Convergence Criterion
    # w is not changed, so we use eVhat$values
    if (isNormal || flag_for_ncj_loop) {
      if (sum(abs((ev_prev - eVhat$values)^2)) <= TolDRes * sum(ev_prev^2)) {
        break
      }

    } else {
      if (sum((w - w_prev)^2) <= TolDRes^2 * sum(w_prev^2)) {
        if (toupper(prior_type) == "CJ") {
          break
        } else {
          flag_for_ncj_loop <- TRUE
        }
      }
    }
  }
  attr(Psihat, "weight.estimated") <- TRUE
  attr(Psihat, "weight") <- w
  attr(Psihat, "noiseCov") <- Vhat

  return(Psihat)

}
