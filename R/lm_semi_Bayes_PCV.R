#' Semiparametric Bayesian shrinkage estimation method for VAR models
#
# lambda can be selected by parameterized cross validation (PCV), and
# lambda_var can be selected by a nonparametric shrinkage type method.
# Additive noise is modeled by multivariate t distribution.
#
#Inputs
#   dof : degree of freedom for multivariate Student t distribution.
#         If Inf, we use multivariate normal distribution.
#         If NULL or a vector of length >=2, we select by K-fold CV.
#   lambda : If NULL or a vector of length >=2,
#         it is selected by PCV (using a default range)
#   lambda_var : If NULL, we estimate it by shrinkage method.
#
# Reference:
#   N. Lee, H. Choi, and S.-H. Kim (2016), CSDA.
#
# 19 Nov. 2017, Namgil Lee, Kangwon National University

lm_semi_Bayes_PCV <- function(Y, X, dof = Inf, lambda = NULL, lambda_var=NULL,
                              prior_type = c("NCJ", "CJ"), num_folds = 5,
                              m0 = ncol(Y)) {

  if (num_folds <= 1)
    stop("Number of folds must be >= 2")

  if (is.null(lambda))
    lambda <- c(1e-3,  seq(1e-2, 1-1e-2, by = 0.01), 1-1e-3, 1-1e-5, 1)

  if (is.null(dof))
    dof <- c(0.2, 0.5, 1, seq(2, 10, by = 2), Inf)

  lenD <- length(dof)
  lenL <- length(lambda)

  d <- ncol(Y)   #size D
  K <- ncol(X)   #size K = D^2p, D(Dp+1)
  p <- K %/% d   #lag order
  N <- nrow(Y)


  #------------- variance -----------------------------------------#
  # Estimation by using (Y, X[,K-d+1 : K])
  #----------------------------------------------------------------#

  # Prepare tsDatc (centered ts data)
  tsDatc <- rbind(X[1:p, (K - d + 1):K], Y)
  tsDatc <- scale(tsDatc, center = TRUE, scale = FALSE)
  lenT <- nrow(tsDatc)  #number of time points

  v1 <- apply(tsDatc, 2, var)  # a row-vector of sample variances

	estimate_lambda_var <- FALSE

	if (is.null(lambda_var)) {
    # Correlation data, w_kjj = (z_kj - z_j)^2, k=1.. lenT, j=1..d
    DatW <- rep(0, d * lenT) 	# a data matrix, transposed
		for (k in 1:lenT) {
			DatW[((k - 1) * d + 1) : (k * d)] <- tsDatc[k, ] ^ 2
		}
		dim(DatW) <- c(d, lenT)

    # Center the correlation data, ie., w_kjj - w_jj
		DatW <- DatW - rowMeans(DatW)	# don't need rep()

    # Compute var_s = var_r + cov_r
		var_r = 1 / ((lenT - 1) ^ 2) * sum(rowSums(DatW ^ 2))
		cov_r = 0
		for (k in 1:(lenT - 1)) {
			for (k2 in (k + 1):lenT) {  # k < k2
				if (k2 - k >= lenT - 1) {
					cov_r <- cov_r + 2 / ((lenT - 1) ^ 2) / lenT *
					  sum( DatW[, -(((lenT - (k2 - k)) + 1):lenT)] * DatW[, -(1:(k2 - k))] )
				} else {
					cov_r <- cov_r + 2 / ((lenT - 1) ^ 2) / lenT *
					  sum( rowSums( DatW[, -(((lenT - (k2 - k)) + 1):lenT)] * DatW[,-(1:(k2 - k))] ) )
				}
			}
		}
		var_s <-  var_r + cov_r

		# Calculate E_vvm2 = sum of squared biases, the denominator
		E_vvm2 = sum( (v1 - median(v1))^2 )

		# Determine lambda_var
		lambda_var = var_s / E_vvm2

		estimate_lambda_var <- TRUE 		#attr(myPsi,"lambda_var.estimated") = TRUE
	}

	# Update variance components
  lambda_var <- max(0, min(1, lambda_var))
  vhat <- (1 - lambda_var) * v1  + lambda_var * median(v1)


	#----------------- correlation ----------------------------------#
	# Estimation using X and Y with each column divided by its std
	#----------------------------------------------------------------#

  # Divide by standard deviation
  Y <- Y / rep( sqrt(v1), each = N)
  X[, (K - d * p + 1):K] <- X[, (K - d * p + 1):K] / rep( sqrt(v1), each = N)

  # Select parameters by K-FOLD CV
  estimate_dof <- FALSE
  estimate_lambda <- FALSE
  if (lenD * lenL > 1) {

    #### Divide sample into num_folds blocks ####
    vidx <- vector("list", num_folds)  # indices for validation data
    tidx <- vector("list", num_folds)  # indices for training data
    idx_s <- sample(N)                  # shupple indices of sample
    bsize <- floor(N / num_folds)        # basic block size
    numAdded <- N - bsize * num_folds  # blocks of size floor(N/K)+1
    numDeflt <- num_folds - numAdded   # blocks of size floor(N/K)
    for (fold in 1:numDeflt) {
      vidx[[fold]] <- idx_s[ (1 + (fold - 1) * bsize) : (fold * bsize) ]
      tidx[[fold]] <- setdiff(idx_s, vidx[[fold]])
    }
    tmpnum = bsize*numDeflt
    for (fold in seq(1, numAdded, length.out = numAdded)) {
      vidx[[fold + numDeflt]] <- idx_s[ (1 + tmpnum + (fold-1)*(bsize+1)) :
                                          (tmpnum + fold*(bsize+1)) ]
      tidx[[fold + numDeflt]] <- setdiff(idx_s, vidx[[fold+numDeflt]])
    }
    ################################

    #### the K-fold CV procedure ############
    sell <- lambda[1]
    seld <- dof[1]
    selMSE <- 1e10
    for (idD in 1:lenD) {
      dof_curr <- dof[idD]

      # Run PCV: use KCV to estimate lambda
      eta_bar <- 0
      #MSE_ave <- 0#not_used
      for (fold in 1:num_folds) {
          XpTrain <- matrix( X[tidx[[fold]], ] ,length(tidx[[fold]]) )
          XfTrain <- matrix( Y[tidx[[fold]], ] ,length(tidx[[fold]]) )
          XpValid <- matrix( X[vidx[[fold]], ] ,length(vidx[[fold]]) )
          XfValid <- matrix( Y[vidx[[fold]], ] ,length(vidx[[fold]]) )

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
                            (d * K))
          eta_bar <- eta_bar + lgtheta1 / num_folds    #log inv eta
      } # end of for (fold)
      etainv <- (d * K) * exp(eta_bar)
      lambda_curr <- etainv / (etainv + N - 1)

      ####### Compute KCV error #######
      MSE_ave <- 0
      for (fold in 1:num_folds) {
          XpTrain <- matrix( X[tidx[[fold]], ], length(tidx[[fold]]) )
          XfTrain <- matrix( Y[tidx[[fold]], ], length(tidx[[fold]]) )
          XpValid <- matrix( X[vidx[[fold]], ], length(vidx[[fold]]) )
          XfValid <- matrix( Y[vidx[[fold]], ], length(vidx[[fold]]) )

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
      #### end of local K-fold CV procedure ####

    }#end for idD in dof

    dof <- seld
    lambda <- sell

    estimate_dof <- (lenD > 1)
    estimate_lambda <- (lenL > 1 )
  }#end if

  lambda <- max(0, min(1, lambda))

  ###########################

  myPsi  <- shrinkVARcoef(Y = Y, X = X, lambda = lambda, dof = dof,
                          prior_type = prior_type, m0 = m0)
  mySigma <- attr(myPsi, "noiseCov")
  myq  <- attr(myPsi, "weight")

  # re-scale myPsi
  myPsi[(K - d * p + 1):K, ] <- (myPsi[(K - d * p + 1):K, ] / sqrt(vhat))
  myPsi <- myPsi * rep(sqrt(vhat), each = K)
  mySigma <- (mySigma * sqrt(vhat)) * rep(sqrt(vhat), each = d)

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
