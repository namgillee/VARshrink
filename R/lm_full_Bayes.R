#' Full Bayesian Shrinkage Estimation Method for Multivariate Regression
#'
#' Estimate regression coefficients and scale matrix for noise by using
#' Gibbs MCMC algorithm. The function assumes 1) multivariate t-distribution
#' for noise as a sampling distribution, and 2) noninformative priors for
#' regression coefficients and scale matrix for noise.
#'
#' Consider the multivariate regression:
#' \deqn{Y = X \Psi + e, \quad e \sim MVT(0, \nu, \Sigma).}
#' \eqn{\Psi} is a M-by-K matrix of regression coefficients and
#' \eqn{\Sigma} is a K-by-K scale matrix for multivariate t-distribution for
#' noise.
#'
#' Sampling distribution for noise e is multivariate t-distribution with
#' degree of freedom dof and scale matrix \eqn{\Sigma: e \sim MVT(0, \nu,
#' \Sigma)}.
#' The priors are 1) the shrinkage prior for regression coefficients
#' \eqn{\Psi}, and 2) the reference prior or inverse-Wishart for scale matrix
#' \eqn{\Sigma}.
#'
#' The function implements Gibbs MCMC algorithm for estimating regression
#' coefficients Psi and scale matrix Sigma.
#'
#' @param Y An N x K matrix of dependent variables.
#' @param X An N x M matrix of regressors.
#' @param dof Degree of freedom for multivariate t-distribution.
#' If \code{dof = Inf} (default), then multivariate normal distribution is
#' applied and weight vector q is not estimated. If \code{dof = NULL} or
#' \code{dof} <= 0, then \code{dof} and q are estimated automatically.
#' If \code{dof} is a positive number, q is estimated.
#' @param prior_type If "NCJ" (default), use shrinkage-reference prior.
#' If "CJ", use shrinkage-inverse Wishart prior.
#' @param burnincycle,mcmccycle Number of burnin cycles is the number of
#' initially generated sample values to drop. Number of MCMC cycles is the
#' number of generated sample values to compute estimates.
#' @param store_mcmc TRUE to return all MCMC chain of estimated parameters
#' across mcmccycle. Default is TRUE.
#' @return A list object with estimated parameters: Psi, Sigma, dof, delta
#' (delta is the reciprocal of lambda), and lambda.
#' Additional components are se.param (standard error of the parameters) and
#' mcmc.param (parameters for whole mcmc chain, stored if store_mcmc is TRUE).
#' @references S. Ni and D. Sun (2005). Bayesian estimates for vector
#' autoregressive models. Journal of Business & Economic Statistics 23(1),
#' 105-117.
#' @importFrom stats rgamma rnorm runif
#' @importFrom utils capture.output

lm_full_Bayes <- function(Y, X, dof = Inf, prior_type = "NCJ",
                          burnincycle = 1000, mcmccycle = 2000,
                          store_mcmc = TRUE) {
  K <- ncol(Y)
  M <- ncol(X)

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
  Sig <- 0.0001 + diag(1 / rgamma(K, shape = 3, scale = 1 / 4), K)  #inv-gamma

  # Gibbs sampler: burn-in cycles & MCMC cycles
  Psi01 <- Sig01 <- dof01 <- mean_delta <- 0
  Psi01SE <- Sig01SE <- dof01SE <- delta01SE <- 0
  if (isTRUE(store_mcmc)) {
    mcmc.param <- vector("list", length = mcmccycle)
  } else {
    mcmc.param <- NULL
  }
  for (i in 1:(burnincycle + mcmccycle)) {
    # Update phi, delta and Sig:

    # (1) Draw phi #######################
    # Draw phi from a multivariate normal distribution.
    # Suppose that phi ~ N_J(mu_Q, V_Q), with V_Q = UDU',
    #	then, U'phi ~ N_J(U' mu_Q, D).
    if (i == 1) {
      eSig <- eigen(Sig)
      invSig <- eSig$vectors %*% (t(eSig$vectors) / pmax(eSig$values, 1e-20))
    }
    eXX <- eigen(t(X) %*% (X * Q))
    mXY <- t(X) %*% (Y * Q) %*% invSig

    if (prior_type == "CJ") {
      # eigenvalues V_Q
      Dv <- 1 / (kronecker(1 / pmax(eSig$values, 1e-14),
                           pmax(eXX$values, 1e-14) + 1 / delta))
      # mean in U'phi ~ N_J(mu0, Dv)
      mu0 <- as.vector((t(eXX$vectors) / (eXX$values + 1 / delta)) %*% mXY %*%
                         eSig$vectors)
    } else {
      # eigenvalues V_Q
      Dv <- 1 / (kronecker(1 / pmax(eSig$values, 1e-14),
                           pmax(eXX$values, 1e-14)) + 1 / delta)
      # mean in U'phi ~ N_J(mu0, Dv)
      mu0 <- Dv * as.vector(t(eXX$vectors) %*% mXY %*% eSig$vectors)
    }
    phi <- mu0 + sqrt(Dv) * rnorm(J, 0, 1)		# U'phi
    dim(phi) <- c(M, K)
    phi <- eXX$vectors %*% (phi %*% t(eSig$vectors))  # U * U'phi

    # (2) Draw delta #####################
    # Draw delta from InvGamma(J/2-1, x/2)
    if (prior_type == "CJ") {
      x <- sum(t((phi %*% eSig$vectors)^2) / pmax(eSig$values, 1e-14))
    } else {
      x <- sum(phi^2)
    }
    if (x < 1e-20)
      delta <- max(x / J, 1e-30)
    else
      delta <- 1 / rgamma(1, shape = J / 2 - 1, scale = 2 / x) # inverse gamma

    # (3) Draw Sig #######################
    Sk <- Y - X %*% phi
    Sk <- t(Sk) %*% (Sk * Q)					# S_k

    if (prior_type == "CJ") {
      # Draw Sig from inverse Wishart
      Sig <- MCMCpack::riwish(K + N, (K + K + 1) * diag(K) + Sk)
      eSig <- eigen(Sig)
      invSig <- eSig$vectors %*% (t(eSig$vectors) / pmax(eSig$values, 1e-20))

    } else {
      # Draw Sig based on Sun & Ni (2004)
      logLambda <- log(pmax(eSig$values, 1e-20))
      SigStar <- (eSig$vectors) %*% diag(logLambda) %*% t(eSig$vectors)

      zij <- rnorm(K * (K + 1) / 2)					# z_ij are standard normals.
      zij <- zij / sqrt(sum(zij^2))				# v_ij are upper-triangular part of V,
      V <- matrix(0, K, K)						#that are normalized standard normals.
      V[upper.tri(V, diag = TRUE)] <- zij			#
      for (j in 2:K) {
        for (k in 1:(j - 1)) {
          V[j, k] <- V[k, j]						# V is a symmetric metric
        }
      }

      W <- SigStar + rnorm(1) * V
      eW <- eigen(W)
      id_sort <- sort.list(eW$values, decreasing = TRUE)
      eW$values <- eW$values[id_sort]
      eW$vectors <- eW$vectors[, id_sort]
      Cstar <- eW$values

      alpha_k <- N / 2 * sum(logLambda - Cstar) +
        1 / 2 * sum((invSig - eW$vectors %*% diag(1 / exp(Cstar)) %*%
                       t(eW$vectors)) * Sk)
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
    }

    # (4) draw Q from gamma ##############
    if (is.null(dof) || !is.infinite(dof)) {
      #If dof is Inf, then it is a multivariate normal distribution,
      #and do not update Q.
      #If dof is not Inf, then it is a multivariate t-distribution, and update Q
      if (is.infinite(w)) {
        Q <- rep(1, N)
      } else {
        x <- (Y - X %*% phi)
        x <- w + 0.5 * rowSums(x * (x %*% invSig))  # x : rate, beta, 1/scale

        #gamma distribution with alpha=(0.5(nu + K)), scale = 1/x
        Q <- rgamma(N, shape = (w + 0.5 * K), scale = 1) / x
      }
    }

    # (5) draw w by MCMC #################
    if (isTRUE(estimate_dof)) {
      f_log <- function(x, N, Q) {
        N * x * log(x) + x * sum(log(Q), na.rm = TRUE) - N * lgamma(x) -
          (1 + sum(Q)) * x
      }
      f_dif <- function(x, N, Q) {
        N * log(x) + N + sum(log(Q), na.rm = TRUE) - N * digamma(x) -
          (1 + sum(Q))
      }
      capture.output({
        w <- ars::ars(n = 1, f = f_log, fprima = f_dif, x = 10, m = 1,
                      lb = TRUE, xlb = 1e-15, N = N, Q = Q)
      })
    }

    # update Psi01, Sig01 #########
    if (i >= (burnincycle + 1)) {

      Psi01 <- Psi01 + phi / mcmccycle
      Sig01 <- Sig01 + Sig / mcmccycle
      mean_delta <- mean_delta + delta / mcmccycle
      dof01 <- dof01 + w * 2 / mcmccycle

      # Note that, SE^2 = S^2/n = sum(x_i^2/n/(n-1)) - xbar^2/(n-1)
      Psi01SE <- Psi01SE + phi^2 / mcmccycle / (mcmccycle - 1)
      Sig01SE <- Sig01SE + Sig^2 / mcmccycle / (mcmccycle - 1)
      delta01SE <- delta01SE + delta^2 / mcmccycle / (mcmccycle - 1)
      if (estimate_dof)
        dof01SE <- dof01SE + (w * 2)^2 / mcmccycle / (mcmccycle - 1)

      # Store parameters across MCMC cycle
      if (isTRUE(store_mcmc)) {
        mcmc.param[[i - burnincycle]] <- list(
          Psi = matrix(phi, M, K),
          Sigma = Sig,
          dof.estimated = estimate_dof,
          dof = if (isTRUE(estimate_dof)) w * 2 else dof,
          delta = delta,
          lambda = 1 / delta,
          lambda.estimated = TRUE
        )
        mcmc_local_dof <- mcmc.param[[i - burnincycle]]$dof
        if (is.infinite(mcmc_local_dof)) {
          mcmc_local_q <- rep(1, N)
        } else {
          x <- (Y - X %*% phi)
          x <- 0.5 * mcmc_local_dof + 0.5 * rowSums(x * t(solve(Sig, t(x))))
          mcmc_local_q <- (0.5 * mcmc_local_dof + 0.5 * K - 1) / pmin(x, 1e10)
        }
        mcmc.param[[i - burnincycle]]$q <- mcmc_local_q
      } # End of if (isTRUE(store_mcmc))
    } # End of if (i >= (burnincycle + 1))
  } # End of for (i in 1:(burnincycle + mcmccycle))
  #######

  # Collect return values
  mean_Psi <- matrix(Psi01, M, K)   # Psi01 matrix
  se.mean_Psi <- matrix(sqrt(Psi01SE - Psi01^2 / (mcmccycle - 1)), M, K)
  mean_Sigma <- Sig01
  se.mean_Sigma <- sqrt(Sig01SE - Sig01^2 / (mcmccycle - 1))
  if (estimate_dof) {
    mean_dof <- dof01
    se.mean_dof <- sqrt(dof01SE - dof01^2 / (mcmccycle - 1))
  } else {
    mean_dof <- dof
    se.mean_dof <- 0
  }
  # Estimate Q values by the mode of posterior, (shape-1)*scale
  # (see, lines 175--179)
  if (is.infinite(mean_dof)) {
    mode_q <- rep(1, N)
    se.mode_q <- rep(0, N)
  } else {
    x <- (Y - X %*% mean_Psi)
    x <- 0.5 * mean_dof + 0.5 * rowSums(x * t(solve(mean_Sigma, t(x))))
    mode_q <- (0.5 * mean_dof + 0.5 * K - 1) / pmin(x, 1e10)
    se.mode_q <- sqrt(0.5 * mean_dof + 0.5 * K) / pmin(x, 1e10)
  }

  res <- list(
    Psi = mean_Psi,
    Sigma = mean_Sigma,
    dof.estimated = estimate_dof,
    dof = mean_dof,
    delta = mean_delta,
    lambda = 1 / mean_delta,
    lambda.estimated = TRUE,
    q = mode_q
  )

  res$se.param <- list(
    Psi = se.mean_Psi,
    Sigma = se.mean_Sigma,
    dof = se.mean_dof,
    delta = sqrt(delta01SE - mean_delta^2 / (mcmccycle - 1)),
    lambda = res$se.param$delta / mean_delta^2,
    q = se.mode_q
  )

  res$mcmc.param <- mcmc.param

  res
}
