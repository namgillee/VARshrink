#' K-fold Cross Validation for Selection of Shrinkage Parameters of
#' Semiparametric Bayesian Shrinkage Estimator for Multivariate Regression
#'
#' Estimate regression coefficients and scale matrix for noise by using
#' semiparametric Bayesian shrinkage estimator, whose shrinkage parameters
#' are selected by k-fold cross validation (KCV).
#'
#' The shrinkage parameters, lambda and lambda_var, for the semiparametric
#' Bayesian shrinkage estimator are selected by KCV. See help(lm_semi_Bayes_PCV)
#' for details about semiparametric Bayesian estimator.
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param dof Degree of freedom for multivariate t-distribution.
#' If dof = Inf (default), then multivariate normal distribution is applied and
#' weight vector q is not estimated. If dof = NULL or a numeric vector,
#' then dof is selected by k-fold CV automatically and q is estimated.
#' @param lambda If NULL or a vector of length >=2, it is selected by KCV.
#' @param lambda_var If NULL or a vector of length >=2, it is selected by KCV.
#' @param prior_type "NCJ" for non-conjugate prior and "CJ" for conjugate
#' prior for scale matrix Sigma.
#' @param num_folds Number of folds for KCV.
#' @param m0 A hyperparameter for inverse Wishart distribution for Sigma
#' @references N. Lee, H. Choi, and S.-H. Kim (2016). Bayes shrinkage
#' estimation for high-dimensional VAR models with scale mixture of normal
#' distributions for noise. Computational Statistics & Data Analysis 101,
#' 250-276. doi: 10.1016/j.csda.2016.03.007
#' @importFrom stats var median
#
# Last modified: 01 Dec. 2017, Namgil Lee @ Kangwon National University

lm_ShVAR_KCV <- function(Y, X, dof = Inf, lambda = NULL, lambda_var = NULL,
                         prior_type = c("NCJ", "CJ"), num_folds = 5,
                         m0 = ncol(Y)) {

  if (num_folds <= 1)
    stop("Number of folds must be >= 2")

  if (is.null(lambda))
    lambda <- c(1e-3,  seq(1e-2, 1 - 1e-2, by = 0.01), 1 - 1e-3, 1 - 1e-5, 1)

  if (is.null(lambda_var))
    lambda_var <- c(1e-3,  seq(1e-2, 1 - 1e-2, by = 0.01), 1 - 1e-3, 1 - 1e-5,
                    1)

  if (is.null(dof))
    dof <- c(0.2, 0.5, 1, seq(2, 10, by = 2), Inf)

  lenLV <- length(lambda_var)
  lenL <- length(lambda)
  lenD <- length(dof)

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

  estimate_lambda_var <- (lenLV >= 2)
  estimate_lambda <- (lenL >= 2)
  estimate_dof <- (lenD >= 2)

  # Extract variance
  #tsDat <- rbind(X[1:p,(M-K+1):M], Y)
  #v1 <- apply(tsDat, 2, var)  # a row-vector of sample variances
  #
  # Divide by standard deviation
  #Note that this standardization process prevents
  #lambda_var to be included in the k-fold CV process,
  #which means that lambda_var cannot be estimated adaptively.
  #Y <- Y / rep( sqrt(v1), each = N)
  #X[,(M-K*p+1):M] <- X[,(M-K*p+1):M] / rep( sqrt(v1), each = N)

  # Select parameters by k-FOLD CV
  if (lenLV * lenL * lenD > 1) {
    #### Divide sample into num_folds blocks ####
    vidx <- vector("list", num_folds)  # indices for validation data
    tidx <- vector("list", num_folds)  # indices for training data
    idx_s <- sample(N)                 # shupple indices of sample
    bsize <- floor(N/num_folds);   # basic block size
    numAdded <- N - bsize * num_folds # blocks of size floor(N/M)+1
    numDeflt <- num_folds - numAdded # blocks of size floor(N/M)
    for (fold in 1:numDeflt) {
      vidx[[fold]] <- idx_s[ (1 + (fold - 1) * bsize):(fold * bsize)]
      tidx[[fold]] <- setdiff(idx_s, vidx[[fold]])
    }
    tmpnum <- bsize * numDeflt;
    for (fold in seq(1,numAdded,length.out=numAdded)) {
      vidx[[fold+numDeflt]] <- idx_s[ (1 + tmpnum + (fold - 1) * (bsize + 1)) :
                                       (tmpnum + fold * (bsize + 1)) ]
      tidx[[fold+numDeflt]] <- setdiff(idx_s, vidx[[fold + numDeflt]])
    }
    ################################

    #### the k-fold CV procedure ############
    sellv <- lambda_var[1]
    sell <- lambda[1]
    seld <- dof[1]
    selMSE <- 1e10
    for (idD in 1:lenD) {
        dof_curr <- dof[idD]

        #how to select lambda_var?
        #For each (lambda, lambda_var), compute MSE_ave, which is
        #the average of k-fold CV prediction errors.
        #Note that Psihat is computed only once for each lambda and fold.
        MSE_ave <- matrix(0, nrow = lenL, ncol = lenLV)
        for (fold in 1:num_folds) {
            XpTrain <- matrix(X[tidx[[fold]], ], length(tidx[[fold]]))
            XfTrain <- matrix(Y[tidx[[fold]], ], length(tidx[[fold]]))
            XpValid <- matrix(X[vidx[[fold]], ], length(vidx[[fold]]))
            XfValid <- matrix(Y[vidx[[fold]], ], length(vidx[[fold]]))

            #### rescale by std
            tsTR <- rbind(XpTrain[1:p, (1 + (p - 1) * K):(p * K)], XfTrain)
            v1TR <- apply(tsTR, 2, var)
            XfTrain <- XfTrain / rep(sqrt(v1TR), each = length(tidx[[fold]]))
            XpTrain[, 1:(p * K)] <- XpTrain[, 1:(p * K)] /
              rep(sqrt(v1TR), each = length(tidx[[fold]]))

            # If dof_curr = Inf, the computation is much easier.
            if (!is.infinite(dof_curr)) {
            		for (idL in 1:lenL) {
            		    lambda_curr <- lambda[idL]
            		    #### estimate Psihat(lambda, dof)
            		    Psihat <- shrinkVARcoef(Y = XfTrain, X = XpTrain,
            		                            lambda = lambda_curr, dof = dof_curr,
            		                            prior_type = prior_type, m0 = m0)
            		    #### estimate lambda_var
            		    for (idLV in 1:lenLV) {
                			lambda_var_curr <- lambda_var[idLV]
                			vhat <- (1 - lambda_var_curr) * v1TR  + lambda_var_curr *
                			  median(v1TR)
                			scaledPsihat <- Psihat
                			scaledPsihat[1:(p * K), ] <-
                			  (scaledPsihat[1:(p * K), ] / sqrt(vhat))
                			scaledPsihat <- scaledPsihat * rep(sqrt(vhat), each = M)

                			pe_k <- sum((XfValid -  XpValid %*% scaledPsihat) ^ 2) /
                			  nrow(XfValid)
                			MSE_ave[idL, idLV] <- MSE_ave[idL, idLV] + pe_k / num_folds
            		    }#end of for(idLV)
            		}#end of for(idL)

            } else {
	              #Estimate Phihat for all theta(lambda) values
                theta <- lambda * (nrow(XfTrain) - 1) / (1 - lambda)

                Xs <- svd(XpTrain)
                Rhs <- t(Xs$u) %*% XfTrain         #r x K

                k <- length(theta)                         #k == lenL
                r <- length(Xs$d)                         #rank
                Div <- Xs$d^2 + rep(theta, each = r * K) #r*K x k
                a <- rep(drop(Xs$d * Rhs), k)/Div
                dim(a) <- c(r, K*k)
                Psihat <- Xs$v %*% a                    #M x K*k

                #### estimate lambda_var
                for (idLV in 1:lenLV) {
                        lambda_var_curr <- lambda_var[idLV]
                        vhat <- (1 - lambda_var_curr) * v1TR  +
                          lambda_var_curr * median(v1TR)
                        scaledPsihat <- Psihat
                        scaledPsihat[1:(p * K),] <-
                          (scaledPsihat[1:(p * K),] / sqrt(vhat))
                        scaledPsihat <- scaledPsihat * rep(sqrt(vhat), each = M)

                        Resid <- rep(XfValid,k) -  XpValid %*% scaledPsihat
                        dim(Resid) <- c(nrow(XpValid) * K, k)
                        pe_k <- colSums(Resid^2) / nrow(XfValid)

                        MSE_ave[, idLV] <- MSE_ave[, idLV] + pe_k / num_folds
                }#end of for(idLV)

            }#end of if(is.infinite)
        }#end of for(fold)

        # Select the (lambda, lambda_var)
        id <- which.min(MSE_ave)
        idLV <- ((id - 1) %/% lenL) + 1  #quotient
        idL <- ((id - 1) %% lenL) + 1     #resid
        mse_curr <- MSE_ave[id]
        lambda_var_curr <- lambda_var[idLV]
        lambda_curr <- lambda[idL]

        if (mse_curr < selMSE) {
            selMSE <- mse_curr
            sellv <- lambda_var_curr
            sell  <- lambda_curr
            seld <- dof_curr
        }
    }#end of for(idD)

    lambda_var <- sellv
    lambda <- sell
    dof <- seld

  }#end of if(lenD*lenL*lenLV)

  lambda <- max(0, min(1, lambda))
  lambda_var <- max(0, min(1, lambda_var))

  ###########################
  # calculate the coefficient matrix
  tsTR <- rbind(X[1:p, (1 + (p - 1) * K):(p * K)], Y)
  v1TR <- apply(tsTR, 2, var)
  Y <- Y / rep(sqrt(v1TR), each = N)
  X[, 1:(p * K)] <- X[, 1:(p * K)] / rep(sqrt(v1TR), each = N)

  myPsi  <- shrinkVARcoef(Y = Y, X = X, lambda = lambda,
	                        dof = dof, prior_type = prior_type, m0 = m0)

  mySigma <- attr(myPsi, "noiseCov")
  myq <- attr(myPsi, "weight")

  # re-scale myPsi
  vhat <- (1 - lambda_var) * v1TR  + lambda_var * median(v1TR)
  myPsi[1:(p * K), ] <- (myPsi[1:(p * K), ] / sqrt(vhat))
  myPsi <- myPsi * rep(sqrt(vhat), each = M)
  mySigma <- (mySigma * sqrt(vhat)) * rep(sqrt(vhat), each = K)

  # Return values
  res <- NULL
  res$Psi <- myPsi
  res$Sigma <- mySigma
  res$dof <- dof
  res$q <- myq
  res$lambda <- lambda
  res$lambda.estimated <- estimate_lambda
  res$lambda_var <- lambda_var
  res$lambda_var.estimated <- estimate_lambda_var
  res$dof.estimated <- estimate_dof

  res
}
