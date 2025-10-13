#' @importFrom stats update quantile
h_boot <- function(x, n.ahead, runs, ortho, cumulative, impulse, response,
                   ci, seed, y.names) {
  if (!(is.null(seed)))
    set.seed(abs(as.integer(seed)))
  if (inherits(x, "varest")) {
    VAR <- eval.parent(x)
  } else if (inherits(x, "svarest")) {
    VAR <- eval.parent(x$var)
  } else {
    stop("Bootstrap not implemented for this class.\n")
  }
  p <- VAR$p
  K <- VAR$K
  obs <- VAR$obs
  total <- VAR$totobs
  B <- Bcoef_sh(VAR)
  BOOT <- vector("list", runs)
  ysampled <- matrix(0, nrow = total, ncol = K)
  colnames(ysampled) <- names(VAR$varresult)
  Zdet <- NULL
  if (ncol(VAR$datamat) > (K * (p + 1))) {
    Zdet <- as.matrix(VAR$datamat[, (K * (p + 1) + 1):ncol(VAR$datamat)])
  }
  resorig <- scale(resid(VAR), scale = FALSE)
  B <- Bcoef_sh(VAR)
  for (i in 1:runs) {
    booted <- sample(c(1:obs), replace = TRUE)
    resid <- resorig[booted, ]
    lasty <- c(t(VAR$y[p:1, ]))
    ysampled[c(1:p), ] <- VAR$y[c(1:p), ]
    for (j in 1:obs) {
      lasty <- lasty[1:(K * p)]
      Z <- c(lasty, Zdet[j, ])
      ysampled[j + p, ] <- B %*% Z + resid[j, ]
      lasty <- c(ysampled[j + p, ], lasty)
    }
    varboot <- update(VAR, y = ysampled)
    if (inherits(x, "svarest")) {
      varboot <- update(x, x = varboot)
    }
    BOOT[[i]] <- h_irf(x = varboot, n.ahead = n.ahead, ortho = ortho,
                       cumulative = cumulative, impulse = impulse,
                       response = response, y.names = y.names)
  }
  lower <- ci / 2
  upper <- 1 - ci / 2
  mat.l <- matrix(NA, nrow = n.ahead + 1, ncol = length(response))
  mat.u <- matrix(NA, nrow = n.ahead + 1, ncol = length(response))
  Lower <- list()
  Upper <- list()
  idx1 <- length(impulse)
  idx2 <- length(response)
  idx3 <- n.ahead + 1
  temp <- rep(NA, runs)
  for (j in 1:idx1) {
    for (m in 1:idx2) {
      for (l in 1:idx3) {
        for (i in 1:runs) {
          if (idx2 > 1) {
            temp[i] <- BOOT[[i]][[j]][l, m]
          } else {
            temp[i] <- matrix(BOOT[[i]][[j]])[l, m]
          }
        }
        mat.l[l, m] <- quantile(temp, lower, na.rm = TRUE)
        mat.u[l, m] <- quantile(temp, upper, na.rm = TRUE)
      }
    }
    colnames(mat.l) <- response
    colnames(mat.u) <- response
    Lower[[j]] <- mat.l
    Upper[[j]] <- mat.u
  }
  names(Lower) <- impulse
  names(Upper) <- impulse
  result <- list(Lower = Lower, Upper = Upper)
  return(result)
}
