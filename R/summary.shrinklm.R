#' Summary method for class "shrinklm"
#'
#' Class "shrinklm" inherits the class "lm", and it extends
#' the "lm" class to incorporate shrinkage estimates with
#' effective number of parameter.
#' @param object An object of class "shrinklm"
#' @param correlation If TRUE, the correlation matrix of the
#' the estimated coefficients is returned and printed.
#' @param symbolic.cor If TRUE, print the correlations in a symbolic form
#' rather than as numbers
#' @importFrom stats coef var pt
#' @export
summary.shrinklm <- function (object, correlation = FALSE,
                              symbolic.cor = FALSE, ...) {
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(z))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <-
      list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(z, "lm"))
    warning("calling summary.lm(<fake-lm-object>) ...")
  if (is.na(z$df.residual))
    warning(paste("residual degrees of freedom in object",
                  "suggest this is not an \"lm\" fit"))
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  n <- length(r)
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept"))
      sum( (f - mean(f))^2 )
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
      1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- (abs(z$svd$d) >= 1e-14)
  R <- z$svd$v[, p1, drop = FALSE] %*% (
    (z$svd$d[p1]^2 / (z$svd$d[p1]^2 + z$lambda0)^2) *
      t(z$svd$v[, p1, drop = FALSE])
  )  ## class 'shrinklm' replaces QR with SVD to compute R ~ (X'X)^{-1} ##
  se <- if (nrow(R) < length(z$coefficients)) {
    # Assume ncol(X) + 1 = length(coef). SE for const is sqrt of
    # var(m_yi + p_i' * m_x) = resvar * (1 + m_x' R m_x)
    m_u <- colMeans(z$svd$u[, p1, drop = FALSE])
    se_const <- sqrt(resvar * (1 + sum(m_u^2 * z$svd$d[p1]^4 /
                                         (z$svd$d[p1]^2 + z$lambda0)^2)))
    c(sqrt(diag(R) * resvar), se_const)
  } else {
    sqrt(diag(R) * resvar)
  }
  est <- z$coefficients
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <-
    cbind(Estimate = est, `Std. Error` = se, `t value` = tval,
          `Pr(>|t|)` = 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  ans$aliased <- is.na(coef(z))
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, p)
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept"))
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
    ans$fstatistic <- c(value = (mss/(n - rdf - df.int)) / resvar,
                        numdf = n - rdf - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <-
    list(dimnames(ans$coefficients)[[1]][1:nrow(R)],
         dimnames(ans$coefficients)[[1]][1:nrow(R)])
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action))
    ans$na.action <- z$na.action
  class(ans) <- "summary.lm"
  ans
}
