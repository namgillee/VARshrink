#' Shrinkage estimation of VAR parameters
#'
#' Shrinkage estimation methods for high-dimensional VAR models.
#' Consider VAR(p) model:
#' \deqn{y_t = A_1 y_{t-1} + \cdots + A_p y_{t-p} + C d_t + e_t,}
#' where \eqn{y_t} is K-dimensional time series,
#' \eqn{d_t} is deterministic regressors, \eqn{e_t} is a noise process, and
#' \eqn{A_1, \ldots, A_p}, and \eqn{C} are coefficient matrices.
#' Exogenous variables can be included additionally as regressors.
#'
#' Shrinkage estimation methods can estimate the coefficients
#' even when the dimensionality K is larger than the
#' number of observations.
#'
#' @param y A T-by-K matrix of endogenous variables
#' @param p Integer for the lag order
#' @param type  Type of deterministic regressors to include.
#' 1) \code{const} - the constant.
#' 2) \code{trend} - the trend.
#' 3) \code{both} - both the constant and the trend.
#' 4) \code{none}  - no deterministic regressors.
#' ***Note: \code{method="ns"} does not accept \code{type="const"} and
#' \code{type="both"} to avoid constant term.
#' @param season An integer value of frequency for inclusion of
#' centered seasonal dummy variables. \code{abs(season)} >= 3.
#' @param exogen A T-by-L matrix of exogenous variables. Default is \code{NULL}.
#' @param method 1) \code{"ridge"} - multivariate ridge regression.
#' 2) \code{"ns"} - a Stein-type nonparametric shrinkage method.
#' 3) \code{"fbayes"} - a full Bayesian shrinkage method using noninformative
#' priors.
#' 4) \code{"sbayes"} - a semiparametric Bayesian shrinkage method using
#' parameterized cross validation.
#' 5) \code{"kcv"} - a semiparametric Bayesian shrinkage method using
#' K-fold cross validation.
#' 6) \code{"ols"} - ordinary least squares method.
#' @param lambda,lambda_var  Shrinkage parameter value(s).
#' Use of this parameter is slightly different for each method:
#' the same value does not imply the same shrinkage estimates.
#' @param dof  Degree of freedom of multivariate t-distribution for noise.
#' Valid only for \code{method = "fbayes"} and \code{method = "sbayes"}.
#' \code{dof = Inf} means multivariate normal distribution.
#' @param ... Extra arguments to pass to a specific function of the
#' estimation method. For example, burnincycle and mcmccycle are for
#' \code{"fbayes"}.
#' @return An object of class "varshrinkest" with the components:
#' varresult, datamat, y, type, p, K, obs,
#' totobs, restrictions, method, lambda, call.
#' The class "varshrinkest" inherits the class "varest" in the package vars.
#' @import vars
#' @examples
#' data(Canada, package = "vars")
#' y <- diff(Canada)
#' VARshrink(y, p = 2, type = "const", method = "ridge")
#' @export
VARshrink  <- function(y, p = 1, type = c("const", "trend", "both", "none"),
                       season = NULL, exogen = NULL,
                       method =
                         c("ridge", "ns", "fbayes", "sbayes", "kcv", "ols"),
                       lambda = NULL, lambda_var = NULL, dof = Inf, ...) {
  cl <- match.call()
  type <- match.arg(type)
  method <- match.arg(method)
  y <- as.matrix(y)
  totobs <- nrow(y)  # Total number of observations
  K <- ncol(y)       # Dimension of output response
  N <- totobs - p    # Sample size

  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if (K < 2)
    stop("The matrix 'y' should contain at least two variables.\n")
  if (is.null(tsnames <- colnames(y))) {
    tsnames <- paste("y", 1:K, sep = "")
    colnames(y) <- tsnames
    warning(paste("No column names supplied in y, using:",
                  paste(tsnames, collapse = ", "), ", instead.\n"))
  }
  if (totobs <= (p + 1))
    stop("Number of total observations must be > p+1\n")
  if (method == "ns") {
    if (identical(type, "const")) {
      type <- "none"
      warning("'ns' method does not allow type='const'.. changed to 'none'.")
    } else if (identical(type, "both")) {
      type <- "trend"
      warning("'ns' method does not allow type='both'.. changed to 'trend'.")
    }
  }

  #========== Build Data Matrices: datX, datY ============

  datY <- y[(p + 1):totobs, ] #N-by-K
  colnames(datY) <- tsnames

  if (identical(type, "const")) {
    M <- K * p + 1
    datX <- matrix(1, N, M)
    for (h in 1:p) {
      datX[, (1 + (h - 1) * K):(h * K)] <- y[(p + 1 - h):(totobs - h), ]
    }
    colnames(datX) <- c(paste(rep(tsnames, times = p), ".l",
                              rep(1:p, each = K), sep = ""), "const")
  } else if (identical(type, "trend")) {
    M <- K * p + 1
    datX <- matrix(1, N, M)
    for (h in 1:p) {
      datX[, (1 + (h - 1) * K):(h * K)] <- y[(p + 1 - h):(totobs - h), ]
    }
    datX[, M] <- seq(p + 1, length.out = N)
    colnames(datX) <- c(paste(rep(tsnames, times = p), ".l",
                              rep(1:p, each = K), sep = ""), "trend")
  } else if (identical(type, "both")) {
    M <- K * p + 2
    datX <- matrix(1, N, M)
    for (h in 1:p) {
      datX[, (1 + (h - 1) * K):(h * K)] <- y[(p + 1 - h):(totobs - h), ]
    }
    datX[, M - 1] <- seq(p + 1, length.out = N)
    colnames(datX) <- c(paste(rep(tsnames, times = p), ".l",
                              rep(1:p, each = K), sep = ""), "trend", "const")
  } else if (identical(type, "none")) {
    M <- K * p
    datX <- matrix(1, N, M) #N-by-M
    for (h in 1:p) {
      datX[, (1 + (h - 1) * K):(h * K)] <- y[(p + 1 - h):(totobs - h), ]
    }
    colnames(datX) <- paste(rep(tsnames, times = p), ".l",
                            rep(1:p, each = K), sep = "")
  } else {
    stop(paste("Unknown type:", type, "\n"))
  }
  if (!(is.null(season)) && (length(season) == 1) && (abs(season) >= 3)) {
    season <- abs(as.integer(season))
    dum <- (diag(season) - 1 / season)[, -season]
    dums <- dum
    while (nrow(dums) < totobs) {
      dums <- rbind(dums, dum)
    }
    dums <- dums[1:totobs, ]
    colnames(dums) <- paste("sd", seq_len(ncol(dums)), sep = "")
    datX <- cbind(datX, dums[-c(1:p), ])
  }
  if (!(is.null(exogen))) {
    exogen <- as.matrix(exogen)
    if (!identical(nrow(exogen), nrow(y))) {
      stop("\nDifferent row size of y and exogen.\n")
    }
    if (is.null(colnames(exogen))) {
      colnames(exogen) <- paste("exo", seq_len(ncol(exogen)), sep = "")
      warning(paste("No column names supplied in exogen, using:",
                    paste(colnames(exogen), collapse = ", "), ", instead.\n"))
    }
    colnames(exogen) <- make.names(colnames(exogen))
    tmp <- colnames(datX)
    datX <- cbind(datX, exogen[-c(1:p), ])
    colnames(datX) <- c(tmp, colnames(exogen))
  }

  #========== Run a Shrinkage Estimation Method ==========

  estim <- list(restrictions = NULL)

  #---------- (1) Ordinary Least Squares ------------
  if (method == "ols") {
    estim <- VAR(y, p = p, type = type, season = season, exogen = exogen, ...)
    estim$lambda <- 0
    estim$lambda.estimated <- FALSE
  }

  #---------- (2) Multivariate Ridge ----------------
  if (method == "ridge") {

    # Compute the coefficient matrix: myPsi (M-by-K)
    # resu_ridge is a list of $Psi, $lambda, $GCV
    resu_ridge  <- lm_multiv_ridge(datY, datX, lambda, ...)
    id_min_gcv <- which.min(resu_ridge$GCV) # Select the minimum GCV
    myPsi <- resu_ridge$Psi[[id_min_gcv]]

    # Update the return value
    estim$varresult <-
      convPsi2varresult(
        Psi = myPsi, Y = datY, X = datX,
        lambda0 = resu_ridge$lambda[id_min_gcv] * N,
        type = type, callstr = cl
      )
    estim$lambda <- resu_ridge$lambda[id_min_gcv]
    estim$lambda.estimated <- as.logical(length(resu_ridge$lambda) > 1)
    estim$GCV <- resu_ridge$GCV[id_min_gcv]

  }
  #---------- (3) Nonparametric Shrinkage -----------
  if (method == "ns") {

    # datY and datX are centered separately by dybar and dxbar.
    dybar <- dxbar <- NULL
    dybar <- colMeans(datY)
    dxbar <- colMeans(datX)
    datY <- datY - rep(dybar, each = N)
    datX <- datX - rep(dxbar, each = N)

    # Estimate covariance matrix S_Z
    Z <- cbind(datX, datY)
    if (is.null(lambda)) {
      if (is.null(lambda_var)) {
        SZ <- corpcor::cov.shrink(Z, verbose = FALSE, ...)
      } else {
        SZ <- corpcor::cov.shrink(Z, lambda.var = lambda_var, verbose = FALSE,
                                  ...)
      }
    } else {
      if (is.null(lambda_var)) {
        SZ <- corpcor::cov.shrink(Z, lambda = lambda, verbose = FALSE, ...)
      } else {
        SZ <- corpcor::cov.shrink(Z, lambda = lambda, lambda.var = lambda_var,
                                  verbose = FALSE, ...)
      }
    }

    # Compute the coefficient matrix Psi by solving linear system
    # Use EVD to avoid singularity
    eigSZ <- eigen(SZ[seq_len(ncol(datX)), seq_len(ncol(datX))])
    idr <- eigSZ$values > 0
    myPsi <- eigSZ$vectors[, idr] %*%
      (1 / eigSZ$values[idr] * t(eigSZ$vectors[, idr])) %*%
      SZ[seq_len(ncol(datX)), (ncol(datX) + 1):ncol(Z)]
    rownames(myPsi) <- colnames(datX)
    colnames(myPsi) <- colnames(datY)

    # Update the return value
    # If type=="const" or "both", then the constant vector will be computed
    # by const' = dybar' - dxbar' %*% Psi and appended to the coefficients.
    estim$varresult <-
      convPsi2varresult(
        Psi = myPsi, Y = datY, X = datX,
        lambda0 = attr(SZ, "lambda") / (1 - attr(SZ, "lambda")) * (N - 1),
        type = type, ybar = dybar, xbar = dxbar,
        callstr = cl
      )
    estim$lambda <- attr(SZ, "lambda")
    estim$lambda.estimated <- attr(SZ, "lambda.estimated")
    estim$lambda_var <- attr(SZ, "lambda.var")
    estim$lambda_var.estimated <- attr(SZ, "lambda.var.estimated")

  }
  #---------- (4) Full Bayesian with Noninformative Priors ----------
  if (method == "fbayes") {

    # Compute the coefficient matrix: myPsi
    # Arguments 'burnincycle' and 'mcmccycle' are included in '...'
    resu_fbayes <- lm_full_Bayes(datY, datX, dof = dof, ...)

    # Update the return value
    estim$lambda <- resu_fbayes$lambda
    estim$lambda.estimated <- resu_fbayes$lambda.estimated
    estim$delta <- resu_fbayes$delta

    estim$Sigma <- resu_fbayes$Sigma
    estim$dof <- resu_fbayes$dof
    estim$dof.estimated <- resu_fbayes$dof.estimated
    estim$mcmc.param <- resu_fbayes$mcmc.param

    myPsi <- resu_fbayes$Psi
    rownames(myPsi) <- colnames(datX)
    colnames(myPsi) <- colnames(datY)
    noise_variances <- if (K >= 2) diag(resu_fbayes$Sigma) else resu_fbayes$Sigma
    estim$varresult <-
      convPsi2varresult(
        Psi = myPsi, Y = datY, X = datX,
        lambda0 = resu_fbayes$lambda,
        type = type,
        Q_values = resu_fbayes$q,
        noise_variances = noise_variances,
        callstr = cl
      )
  }
  #---------- (5) Semi-parametric Bayesian with lambda by P-CV ----------
  if (method == "sbayes") {

    resu_sbayes <-
      lm_semi_Bayes_PCV(datY, datX, dof = dof, lambda = lambda,
                        lambda_var = lambda_var, ...)

    # Update the return value
    estim$lambda <- resu_sbayes$lambda
    estim$lambda.estimated <- resu_sbayes$lambda.estimated
    estim$lambda_var <- resu_sbayes$lambda_var
    estim$lambda_var.estimated <- resu_sbayes$lambda_var.estimated

    estim$Sigma <- resu_sbayes$Sigma
    estim$dof <- resu_sbayes$dof
    estim$dof.estimated <- resu_sbayes$dof.estimated

    myPsi <- resu_sbayes$Psi
    rownames(myPsi) <- colnames(datX)
    colnames(myPsi) <- colnames(datY)
    noise_variances <- if (K >= 2) diag(resu_sbayes$Sigma) else resu_sbayes$Sigma
    estim$varresult <-
      convPsi2varresult(
        Psi = myPsi, Y = datY, X = datX,
        lambda0 = resu_sbayes$lambda / (1 - resu_sbayes$lambda) * (N - 1),
        type = type,
        Q_values = resu_sbayes$q,
        noise_variances = noise_variances,
        callstr = cl
      )
  }
  #---------- (6) Semi-parametric Bayesian with lambda by K-CV ----------
  if (method == "kcv") {

    resu_kcv <- lm_ShVAR_KCV(datY, datX, dof = dof,
                             lambda = lambda, lambda_var = lambda_var, ...)

    # Update the return value
    estim$lambda <- resu_kcv$lambda
    estim$lambda.estimated <- resu_kcv$lambda.estimated
    estim$lambda_var <- resu_kcv$lambda_var
    estim$lambda_var.estimated <- resu_kcv$lambda_var.estimated

    estim$Sigma <- resu_kcv$Sigma
    estim$dof <- resu_kcv$dof
    estim$dof.estimated <- resu_kcv$dof.estimated

    myPsi <- resu_kcv$Psi
    rownames(myPsi) <- colnames(datX)
    colnames(myPsi) <- colnames(datY)
    noise_variances <- if (K >= 2) diag(resu_kcv$Sigma) else resu_kcv$Sigma
    estim$varresult <-
      convPsi2varresult(
        Psi = myPsi, Y = datY, X = datX,
        lambda0 = resu_kcv$lambda / (1 - resu_kcv$lambda) * (N - 1),
        type = type,
        Q_values = resu_kcv$q,
        noise_variances = noise_variances,
        callstr = cl
      )
  }

  #========== Check return value =========================
  # Components for 'varest'
  if (!with(estim, exists("varresult"))) {
    warning("VAR parameters were not estimated. Check the method.")
  }
  if (is.null(estim$datamat))
    estim$datamat <- data.frame(cbind(datY, datX))
  estim$y        <- y
  estim$type     <- type
  estim$p        <- p
  estim$K        <- K
  estim$obs      <- N
  estim$totobs   <- totobs
  estim$call     <- cl

  # Components for 'varshrinkest'
  estim$method   <- method

  class(estim) <- c("varshrinkest", "varest")

  return(estim)
}
