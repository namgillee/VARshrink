## R code for the article
## "VARshrink: Shrinkage Estimation Methods for Vector Autoregressive Models"
## by Namgil Lee and Sung-Ho Kim
## 6 June 2025

## ----setup, include = FALSE----------------------------------------------
library(VARshrink)


## ---- results = "hide", message=FALSE------------------------------------
set.seed(1000)
myCoef <- list(A = list(matrix(c(0.5, 0, 0, 0.5), 2, 2)), c = c(0.2, 0.7))
myModel <- list(Coef = myCoef, Sigma = diag(0.1^2, 2), dof = Inf)
Y <- simVARmodel(numT = 100, model = myModel, burnin = 10)


## ------------------------------------------------------------------------
resu_estim <- list()
resu_estim$`Ridge regression` <-
  VARshrink(Y, p = 1, type = "const", method = "ridge", lambda = NULL)
resu_estim$`Ridge regression`


## ------------------------------------------------------------------------
summary(resu_estim$`Ridge regression`)


## ------------------------------------------------------------------------
c(AIC = AIC(resu_estim$`Ridge regression`),
  BIC = BIC(resu_estim$`Ridge regression`))


## ------------------------------------------------------------------------
resu_estim$`Nonparametric shrinkage` <-
  VARshrink(Y, p = 1, type = "none", method = "ns", lambda = NULL,
            lambda_var = NULL)
resu_estim$`Nonparametric shrinkage`


## ------------------------------------------------------------------------
summary(resu_estim$`Nonparametric shrinkage`)


## ------------------------------------------------------------------------
resu_estim$`Full Bayes (fixed dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "fbayes", dof = 6,
            burnincycle = 1000, mcmccycle = 2000)
resu_estim$`Full Bayes (fixed dof)`


## ------------------------------------------------------------------------
summary(resu_estim$`Full Bayes (fixed dof)`)


## ------------------------------------------------------------------------
resu_estim$`Full Bayes (estim dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "fbayes", dof = NULL,
            burnincycle = 1000, mcmccycle = 2000)
resu_estim$`Full Bayes (estim dof)`


## ------------------------------------------------------------------------
resu_estim$`Semi Bayes (fixed dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "sbayes", dof = 6, lambda = NULL,
            lambda_var = NULL, prior_type = "NCJ", num_folds = 5, m0 = ncol(Y))
resu_estim$`Semi Bayes (fixed dof)`
summary(resu_estim$`Semi Bayes (fixed dof)`)


## ------------------------------------------------------------------------
resu_estim$`Semi Bayes (estim dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "sbayes", dof = NULL,
            lambda = NULL, lambda_var = NULL, prior_type = "NCJ", num_folds = 5,
            m0 = ncol(Y))
resu_estim$`Semi Bayes (estim dof)`


## ------------------------------------------------------------------------
resu_estim$`K-fold CV (fixed dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "kcv", dof = 6, lambda = NULL,
            lambda_var = NULL, prior_type = "NCJ", num_folds = 5, m0 = ncol(Y))
resu_estim$`K-fold CV (fixed dof)`


## ---- results = "hide"---------------------------------------------------
resu_estim$`K-fold CV (estim dof)` <-
  VARshrink(Y, p = 1, type = "const", method = "kcv", dof = NULL, lambda = NULL,
            lambda_var = NULL, prior_type = "NCJ", num_folds = 5, m0 = ncol(Y))


## ------------------------------------------------------------------------
resu_sse <-
  data.frame(SSE = sapply(resu_estim,
                          function(x) calcSSE_Acoef(Acoef_sh(x), myCoef$A)))


## ---- echo = FALSE-------------------------------------------------------
knitr::kable(round(resu_sse, 3),
             caption = paste0("Sum of squared errors of VAR coefficients ",
                              "estimated by the shrinkage methods."))


## ----diffCanada, fig.cap = paste0("Benchmark dataset obtained by ",
##     "differencing the Canada time series from the \pkg{vars} package")----
data(Canada, package = "vars")
Y <- diff(Canada)
plot(Y, cex.lab = 1.3)


## ------------------------------------------------------------------------
set.seed(1000)
resu_model <-
  array(NA,
        dim = c(5, 2, 3),
        dimnames = list(c("Ridge regression", "Nonparametric shrinkage",
                          "Full Bayes", "Semi Bayes", "K-fold CV"),
                        c("AIC", "BIC"),
                        c("p=1", "p=2", "p=3")))
for (p in 1:3) {
  EstimRidge <- VARshrink(Y, p = p, type = "const", method = "ridge")
  resu_model["Ridge regression", , p] <- c(AIC(EstimRidge), BIC(EstimRidge))

  EstimNS <- VARshrink(Y, p = p, type = "none", method = "ns")
  resu_model["Nonparametric shrinkage", , p] <-
    c(AIC(EstimNS), BIC(EstimNS))

  EstimFB <- VARshrink(Y, p = p, type = "const", method = "fbayes", dof = NULL)
  resu_model["Full Bayes", , p] <- c(AIC(EstimFB), BIC(EstimFB))

  EstimSB <- VARshrink(Y, p = p, type = "const", method = "sbayes",
                       dof = NULL, prior_type = "NCJ")
  resu_model["Semi Bayes", , p] <- c(AIC(EstimSB), BIC(EstimSB))

  EstimKCV <- VARshrink(Y, p = p, type = "const", method = "kcv",
                        dof = NULL, prior_type = "NCJ")
  resu_model["K-fold CV", , p] <- c(AIC(EstimKCV), BIC(EstimKCV))
}


## ----modelcomp, echo = FALSE---------------------------------------------
knitr::kable(round(resu_model, 1),
             caption = paste0("\\label{tab:modelcomp}Information criteria ",
                              "(AIC and BIC) for model comparison."))


## ----pred, fig.cap=paste0("10-step-ahead forecast of the differenced Canada ",
##     "time series using the VAR model estimated by the NS method with ",
##     "lag order $p=2$.")----
plot(predict(VARshrink(Y, p = 2, type = "none", method = "ns")), names = "U")


## ----comparisons_K20, eval=FALSE-----------------------------------------
set.seed(1000)
p <- 1
K <- 20
NumTimePTS <- c(20, 40, 80, 160)
const_vector <- c(rep(0.2, 5), rep(0.7, 15))
Sig <- diag(0.5, K) + matrix(0.5, K, K)
resu_SSE <-
  array(0, dim = c(50, 4, length(NumTimePTS)),
        dimnames = list(1:50, c("Ridge", "NS", "FBayes", "SBayes"), NumTimePTS))
for (idT in seq_along(NumTimePTS)) {
  numT <- NumTimePTS[idT]
  for (r in 1:50) {
    Ad <- createVARCoefs_ltriangular(p = p, K = K, diag_val = 0.6,
                                     num_nonzero = K,
                                     const_vector = const_vector, range_max = 1)
    Md <- list(Coef = Ad, Sigma = Sig, dof = Inf)
    Y <- simVARmodel(numT = numT, model = Md, burnin = 20)
    EstimRG <- VARshrink(Y, p, "const", method = "ridge")
    EstimNS <- VARshrink(Y, p, "none", method = "ns")
    EstimFB <- VARshrink(Y, p, "const", method = "fbayes", dof = NULL)
    EstimSB <- VARshrink(Y, p, "const", method = "sbayes", prior_type = "NCJ")
    resu_SSE[r, 1, idT] <- calcSSE_Acoef(Acoef_sh(EstimRG), Ad$A)
    resu_SSE[r, 2, idT] <- calcSSE_Acoef(Acoef_sh(EstimNS), Ad$A)
    resu_SSE[r, 3, idT] <- calcSSE_Acoef(Acoef_sh(EstimFB), Ad$A)
    resu_SSE[r, 4, idT] <- calcSSE_Acoef(Acoef_sh(EstimSB), Ad$A)
  }
}
op <- par(mfrow = c(2, 2), mar = c(2.2, 4.2, 2.5, 1))
for (idT in seq_along(NumTimePTS)) {
  boxplot(resu_SSE[, , idT], log = "y", cex.lab = 1.5, cex.main = 1.5,
          ylab = "SSE", xlab = "", ylim = c(min(resu_SSE), max(resu_SSE)),
          main = paste("T =", NumTimePTS[idT]))
  grid()
}
par(op)
