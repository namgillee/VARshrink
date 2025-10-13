#' Semiparametric Bayesian Shrinkage Estimation Method for Multivariate
#' Regression
#'
#' Estimate regression coefficients and scale matrix for noise by using
#' a parameterized cross validation (PCV). The function assumes
#' 1) multivariate t-distribution for noise as a sampling distribution,
#' and 2) informative priors for regression coefficients and scale matrix
#' for noise.
#'
#' Consider the multivariate regression:
#' \deqn{Y = X Psi + e, \quad e ~ mvt(0, dof, Sigma).}
#' Psi is a M-by-K matrix of regression coefficients and
#' Sigma is a K-by-K scale matrix for multivariate t-distribution for noise.
#'
#' Sampling distribution for noise e is the multivariate t-distribution with
#' degree of freedom dof and scale matrix Sigma: e ~ mvt(0, dof, Sigma).
#' The priors are informative priors: 1) a shrinkage prior for regression
#' coefficients Psi, and 2) inverse Wishart prior for scale matrix Sigma,
#' which can be either non-conjugate ("NCJ") or conjugate ("CJ") to the
#' shrinkage prior for coefficients Psi.
#'
#' The function implements parameterized cross validation (PCV) for
#' selecting a shrinkage parameter lambda for estimating regression
#' coefficients (0 < lambda <= 1).
#' In addition, the function uses a Stein-type shrinkage method for selecting
#' a shrinkage parameter lambda_var for estimating variances of
#' time series variables.
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param dof Degree of freedom for multivariate t-distribution.
#' If dof = Inf (default), then multivariate normal distribution is applied and
#' weight vector q is not estimated. If dof = NULL or a numeric vector,
#' then dof is selected by k-fold CV automatically and q is estimated.
#' @param lambda If NULL or a vector of length >=2, it is selected by PCV.
#' @param lambda_var If NULL, it is selected by a Stein-type shrinkage method.
#' @param prior_type "NCJ" for non-conjugate prior and "CJ" for conjugate
#' prior for scale matrix Sigma.
#' @param num_folds Number of folds for PCV.
#' @param m0 A hyperparameter for inverse Wishart distribution for Sigma
#' @references N. Lee, H. Choi, and S.-H. Kim (2016). Bayes shrinkage
#' estimation for high-dimensional VAR models with scale mixture of normal
#' distributions for noise. Computational Statistics & Data Analysis 101,
#' 250-276. doi: 10.1016/j.csda.2016.03.007
#' @importFrom stats var median
#
# Last modified: 13 Oct 2025, Namgil Lee @ Kangwon National University

lm_semi_Bayes_PCV <- function(Y, X, dof = Inf, lambda = NULL, lambda_var = NULL,
                              prior_type = c("NCJ", "CJ"), num_folds = 5,
                              m0 = ncol(Y)) {

  if (num_folds <= 1)
    stop("Number of folds must be >= 2")

  if (is.null(lambda))
    lambda <- c(1e-3,  seq(1e-2, 1 - 1e-2, by = 0.01), 1 - 1e-3, 1 - 1e-5, 1)

  if (is.null(dof))
    dof <- c(0.2, 0.5, 1, seq(2, 10, by = 2), Inf)

  lenD <- length(dof)
  lenL <- length(lambda)

  K <- ncol(Y)
  M <- ncol(X)
  N <- nrow(Y)

  # Determine lag order
  for (p in 1:N) {
    col_start <- 1 + p * K
    col_end <- (p + 1) * K
    if (col_end > M) {
      break
    } else if (sum(abs(X[p + 1, col_start:col_end] - X[1, 1:K])) >
                 sqrt(.Machine$double.eps)) {
      break
    }
  }

  #------------- variance -----------------------------------------#
  # Estimation by using centered time series data
  #----------------------------------------------------------------#

  # Prepare tsDatc (centered ts data)
  tsDatc <- rbind(X[1:p, (1 + (p - 1) * K):(p * K)], Y)
  tsDatc <- scale(tsDatc, center = TRUE, scale = FALSE)
  lenT <- nrow(tsDatc)  #number of time points

  v1 <- apply(tsDatc, 2, var)  # a row-vector of sample variances

	estimate_lambda_var <- FALSE

	if (is.null(lambda_var)) {
    # Correlation data, w_kjj = (z_kj - z_j)^2, k=1.. lenT, j=1..K
    DatW <- rep(0, K * lenT) 	# a data matrix, transposed
		for (k in 1:lenT) {
			DatW[((k - 1) * K + 1) : (k * K)] <- tsDatc[k, ] ^ 2
		}
		dim(DatW) <- c(K, lenT)

    # Center the correlation data, ie., w_kjj - w_jj
		DatW <- DatW - rowMeans(DatW)	# don't need rep()

    # Compute var_s = var_r + cov_r
		var_r <- 1 / ((lenT - 1) ^ 2) * sum(rowSums(DatW ^ 2))
		cov_r <- 0
		for (k in 1:(lenT - 1)) {
			for (k2 in (k + 1):lenT) {
				if (k2 - k >= lenT - 1) {
					cov_r <- cov_r + 2 / ((lenT - 1) ^ 2) / lenT *
					  sum(DatW[, -(((lenT - (k2 - k)) + 1):lenT)] * DatW[, -(1:(k2 - k))])
				} else {
					cov_r <- cov_r + 2 / ((lenT - 1) ^ 2) / lenT *
					  sum(rowSums(DatW[, -(((lenT - (k2 - k)) + 1):lenT)] *
					                DatW[, -(1:(k2 - k))]))
				}
			}
		}
		var_s <-  var_r + cov_r

		# Calculate E_vvm2 = sum of squared biases, the denominator
		E_vvm2 <- sum((v1 - median(v1))^2)

		# Determine lambda_var
		lambda_var <- var_s / E_vvm2

		estimate_lambda_var <- TRUE
	}

	# Update variance components
  lambda_var <- max(0, min(1, lambda_var))
  vhat <- (1 - lambda_var) * v1  + lambda_var * median(v1)


	#----------------- correlation ----------------------------------#
	# Estimation using X and Y with each column divided by its std
	#----------------------------------------------------------------#

  # Divide by standard deviation
  Y <- Y / rep(sqrt(v1), each = N)
  X[, 1:(p * K)] <- X[, 1:(p * K)] / rep(sqrt(v1), each = N)

  # Select parameters by k-FOLD CV
  estimate_dof <- FALSE
  estimate_lambda <- FALSE
  if (lenD * lenL > 1) {

    #### Divide sample into num_folds blocks ####
    vidx <- vector("list", num_folds)  # indices for validation data
    tidx <- vector("list", num_folds)  # indices for training data
    idx_s <- sample(N)                  # shuffle indices of sample
    bsize <- floor(N / num_folds)        # basic block size
    numAdded <- N - bsize * num_folds  # blocks of size floor(N/M)+1
    numDeflt <- num_folds - numAdded   # blocks of size floor(N/M)
    for (fold in 1:numDeflt) {
      vidx[[fold]] <- idx_s[(1 + (fold - 1) * bsize) : (fold * bsize)]
      tidx[[fold]] <- setdiff(idx_s, vidx[[fold]])
    }
    tmpnum <- bsize * numDeflt
    for (fold in seq(1, numAdded, length.out = numAdded)) {
      vidx[[fold + numDeflt]] <- idx_s[(1 + tmpnum + (fold - 1) * (bsize + 1)) :
                                         (tmpnum + fold * (bsize + 1))]
      tidx[[fold + numDeflt]] <- setdiff(idx_s, vidx[[fold + numDeflt]])
    }
    ################################

    #### the k-fold CV procedure ############
    sell <- lambda[1]
    seld <- dof[1]
    selMSE <- 1e10
    for (idD in 1:lenD) {
      dof_curr <- dof[idD]

      # Run PCV: use KCV to estimate lambda
      eta_bar <- 0
      for (fold in 1:num_folds) {
        XpTrain <- matrix(X[tidx[[fold]], ], length(tidx[[fold]]))
        XfTrain <- matrix(Y[tidx[[fold]], ], length(tidx[[fold]]))
        XpValid <- matrix(X[vidx[[fold]], ], length(vidx[[fold]]))
        XfValid <- matrix(Y[vidx[[fold]], ], length(vidx[[fold]]))

        # train and test : select lambda^*_fold
        lambd1 <- 1
        mse1 <- 1e10
        for (idL in 1:lenL) {
          lambd2 <- lambda[idL]
          Psihat <- shrinkVARcoef(Y = XfTrain, X = XpTrain,
                                  lambda = lambd2, dof = dof_curr,
                                  prior_type = prior_type, m0 = m0)
          pe_k <- sum((XfValid -  XpValid %*% Psihat)^2) / nrow(XfValid)

          #Update minimum PE
          if (pe_k < mse1) {
            lambd1 <- lambd2
            mse1 <- pe_k
          }
        }
        lgtheta1 <- log(lambd1 * (nrow(XpTrain) - 1) / (1 - lambd1) /
                          (K * M))
        eta_bar <- eta_bar + lgtheta1 / num_folds    #log inv eta
      } # end of for (fold)
      etainv <- (K * M) * exp(eta_bar)
      lambda_curr <- etainv / (etainv + N - 1)

      ####### Compute KCV error #######
      MSE_ave <- 0
      for (fold in 1:num_folds) {
        XpTrain <- matrix(X[tidx[[fold]], ], length(tidx[[fold]]))
        XfTrain <- matrix(Y[tidx[[fold]], ], length(tidx[[fold]]))
        XpValid <- matrix(X[vidx[[fold]], ], length(vidx[[fold]]))
        XfValid <- matrix(Y[vidx[[fold]], ], length(vidx[[fold]]))

        Psihat <- shrinkVARcoef(Y = XfTrain, X = XpTrain,
                                lambda = lambda_curr, dof = dof_curr,
                                prior_type = prior_type, m0 = m0)
        pe_fold <- sum((XfValid -  XpValid %*% Psihat)^2) / nrow(XfValid)
        MSE_ave <- MSE_ave + pe_fold / num_folds
      }
      ############################

      if (MSE_ave < selMSE) {
        #Update seld, sell
        selMSE <- MSE_ave
        sell <- lambda_curr
        seld <- dof_curr
      }
      #### end of local k-fold CV procedure ####

    }#end for idD in dof

    dof <- seld
    lambda <- sell

    estimate_dof <- (lenD > 1)
    estimate_lambda <- (lenL > 1)
  }#end if

  lambda <- max(0, min(1, lambda))

  ###########################

  myPsi  <- shrinkVARcoef(Y = Y, X = X, lambda = lambda, dof = dof,
                          prior_type = prior_type, m0 = m0)
  mySigma <- attr(myPsi, "noiseCov")
  myq  <- attr(myPsi, "weight")

  # re-scale myPsi
  myPsi[1:(p * K), ] <- (myPsi[1:(p * K), ] / sqrt(vhat))
  myPsi <- myPsi * rep(sqrt(vhat), each = M)
  mySigma <- (mySigma * sqrt(vhat)) * rep(sqrt(vhat), each = K)

  # Return values
  res <- NULL
  res$Psi <- myPsi
  res$Sigma <- mySigma
  res$dof <- dof
  res$lambda <- lambda
  res$lambda.estimated <- estimate_lambda
  res$lambda_var <- lambda_var
  res$lambda_var.estimated <- estimate_lambda_var
  res$dof.estimated <- estimate_dof
  res$q <- myq

  res
}
