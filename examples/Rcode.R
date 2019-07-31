#### R code for whole examples in the paper ####

#--- Section 3. Shrinkage Estimation Methods ---#
# Simulate VAR processes with multivariate normal distribution,
# and run Shrinkage Estimation methods.
#
# # Source all files from directories
# file.sources = list.files(c("./R"),
#                           pattern="*.R$", full.names=TRUE,
#                           ignore.case=TRUE)
# sapply(file.sources,source,.GlobalEnv)


#install("../VARshrink")
#library(VARshrink)
library(devtools)
devtools::load_all(".") #Ctrl + Shift + L

# Random number generator
set.seed(100)

# Create a VAR model with p=1 and multivariate normal distribution (dof = Inf)
myCoef = list(A = list(matrix(c(0.5, 0, 0, 0.5), 2, 2)), c = c(0.2, 0.7))
myModel = VARparam(Coef = myCoef, Sigma = diag(0.1^2, 2), dof = Inf)

# Simulate VAR processes with multivariate normal distribution
Y = simVARmodel(numT = 100, model = myModel, burnin = 10)

# Run shrinkage estimation methods
# 1) Multivariate Ridge Regression
EstimRidge = VARshrink(Y, p = 1, type = 'const', method = 'ridge',
                       lambda = NULL)
EstimRidge
summary(EstimRidge)

# 2) A Nonparametric Shrinkage Method by Opgen-Rein & Strimmer (2007)
EstimNS = VARshrink(Y, p = 1, type = 'const', method = 'ns',
                    lambda = NULL, lambda_var = NULL)
EstimNS
summary(EstimNS)

# 3-1) A Full Bayesian Method Using Gibbs MCMC (known fixed dof) by Ni & Sun (2005)
EstimFB1 = VARshrink(Y, p = 1, type = 'const', method = 'fbayes', dof = 6,
                     burnincycle = 1000, mcmccycle = 2000)
EstimFB1
summary(EstimFB1)

# 3-2) (unknown dof)
EstimFB2 = VARshrink(Y, p = 1, type = 'const', method = 'fbayes', dof = NULL,
                     burnincycle = 1000, mcmccycle = 2000)
EstimFB2
summary(EstimFB2)

# 4-1) A semiparametric Shrinkage method using PCV (known fixed dof) by Lee, Choi & Kim (2016)
EstimSB1 = VARshrink(Y, p = 1, type = 'const', method = 'sbayes', dof = 6,
                     lambda = NULL, lambda_var = NULL, prior_type = 'NCJ',
                     num_folds = 5, m0 = ncol(Y))
EstimSB1
summary(EstimSB1)

# 4-2) (unknown dof)
EstimSB2 = VARshrink(Y, p = 1, type = 'const', method = 'sbayes', dof = NULL,
                    lambda = NULL, lambda_var = NULL, prior_type = 'NCJ',
                    num_folds = 5, m0 = ncol(Y))
EstimSB2
summary(EstimSB2)

# 5) A semiparametric shrinkage method using K-fold Cross Validation
# 5-1) (known fixed dof)
EstimKCV1 = VARshrink(Y, p = 1, type = 'const', method = 'kcv', dof = 6,
                      lambda = NULL, lambda_var = NULL, prior_type = 'NCJ',
                      num_folds = 5, m0 = ncol(Y))
EstimKCV1
summary(EstimKCV1)

# 5-2) (unknown dof)
EstimKCV2 = VARshrink(Y, p = 1, type = 'const', method = 'kcv', dof = NULL,
                 lambda = NULL, lambda_var = NULL, prior_type = 'NCJ',
                 num_folds = 5, m0 = ncol(Y))
EstimKCV2
summary(EstimKCV2)

print('---- VARshrink: Shrinkage Estimation for VAR models ----')
print('  Results of sum of squared errors of estimators     ')
print('')
print(sprintf('Ridge regression       : %.4f', calcSSE_VARcoef(myModel, EstimRidge$varparam)))
print(sprintf('Nonparametric shrinkage: %.4f', calcSSE_VARcoef(myModel, EstimNS$varparam)))
print(sprintf('Full Bayes(fixed dof)  : %.4f', calcSSE_VARcoef(myModel, EstimFB1$varparam)))
print(sprintf('Full Bayes(estim dof)  : %.4f', calcSSE_VARcoef(myModel, EstimFB2$varparam)))
print(sprintf('Semi Bayes(fixed dof)  : %.4f', calcSSE_VARcoef(myModel, EstimSB1$varparam)))
print(sprintf('Semi Bayes(estim dof)  : %.4f', calcSSE_VARcoef(myModel, EstimSB2$varparam)))
print(sprintf('K-fold CV (fixed dof)  : %.4f', calcSSE_VARcoef(myModel, EstimKCV1$varparam)))
print(sprintf('K-fold CV (estim dof)  : %.4f', calcSSE_VARcoef(myModel, EstimKCV2$varparam)))
print('-----------------------------------------------------')


#--- Section 4.1 Benchmark Data (Example_02.R) ---#
#Canada data
data(Canada, package = "vars")
Y = diff(Canada)
plot(Y, cex.lab = 1.3)
set.seed(100)

#type = 'const', method = 'ridge'
EstimRidge10 = VARshrink(Y, p = 1, type = 'const', method = 'ridge')
EstimRidge20 = VARshrink(Y, p = 2, type = 'const', method = 'ridge')
EstimRidge30 = VARshrink(Y, p = 3, type = 'const', method = 'ridge')
#type = 'const', method = 'ns'
EstimNS10 = VARshrink(Y, p = 1, type = 'const', method = 'ns')
EstimNS20 = VARshrink(Y, p = 2, type = 'const', method = 'ns')
EstimNS30 = VARshrink(Y, p = 3, type = 'const', method = 'ns')
#type = 'const', method = 'fbayes', dof = NULL
EstimFB10 = VARshrink(Y, p = 1, type = 'const', method = 'fbayes',
                      dof = NULL, burnincycle = 1000, mcmccycle = 2000)
EstimFB20 = VARshrink(Y, p = 2, type = 'const', method = 'fbayes',
                      dof = NULL, burnincycle = 1000, mcmccycle = 2000)
EstimFB30 = VARshrink(Y, p = 3, type = 'const', method = 'fbayes',
                      dof = NULL, burnincycle = 1000, mcmccycle = 2000)
#type = 'const', method = 'sbayes', dof = c(2,5,Inf)
EstimSB10 = VARshrink(Y, p = 1, type = 'const', method = 'sbayes',
                      dof = c(2,5,Inf), lambda = NULL, lambda_var = NULL,
                      prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimSB20 = VARshrink(Y, p = 2, type = 'const', method = 'sbayes',
                      dof = c(2,5,Inf), lambda = NULL, lambda_var = NULL,
                      prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimSB30 = VARshrink(Y, p = 3, type = 'const', method = 'sbayes',
                      dof = c(2,5,Inf), lambda = NULL, lambda_var = NULL,
                      prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))

#Additionally, the semiparametric Bayes method using K-fold CV
EstimKCV20 = VARshrink(Y, p = 2, type = 'const', method = 'kcv', dof = c(2,5,Inf),
                      lambda = NULL, lambda_var = NULL, prior_type = 'NCJ',
                      num_folds = 5, m0 = ncol(Y))

cat("======== Model Comparison ========\n")
cat("AIC:\n")
cat("Method:  p = 1    p = 2    p = 3  \n")
cat("ridge : ", round(AIC(EstimRidge10),1), "  ", round(AIC(EstimRidge20),1),
    "  ", round(AIC(EstimRidge30),1), "\n")
cat("ns    : ", round(AIC(EstimNS10),1), "  ", round(AIC(EstimNS20),1),
    "  ", round(AIC(EstimNS30),1), "\n")
cat("fbayes: ", round(AIC(EstimFB10),1), "  ", round(AIC(EstimFB20),1),
    "  ", round(AIC(EstimFB30),1), "\n")
cat("sbayes: ", round(AIC(EstimSB10),1), "  ", round(AIC(EstimSB20),1),
    "  ", round(AIC(EstimSB30),1), "\n")
cat("BIC:\n")
cat("Method:  p = 1    p = 2    p = 3  \n")
cat("ridge : ", round(BIC(EstimRidge10),1), "  ", round(BIC(EstimRidge20),1),
    "  ", round(BIC(EstimRidge30),1), "\n")
cat("ns    : ", round(BIC(EstimNS10),1), "  ", round(BIC(EstimNS20),1),
    "  ", round(BIC(EstimNS30),1), "\n")
cat("fbayes: ", round(BIC(EstimFB10),1), "  ", round(BIC(EstimFB20),1),
    "  ", round(BIC(EstimFB30),1), "\n")
cat("sbayes: ", round(BIC(EstimSB10),1), "  ", round(BIC(EstimSB20),1),
    "  ", round(BIC(EstimSB30),1), "\n")

# Best estimator of VAR parameters:
EstimNS30


#--- Section 4.2 Comparative Simulations (comparisons_d20.R) ---#
set.seed(100)
methods <- c("ridge", "ns", "fbayes", "sbayes")
p <- 1; d <- 20; NumTimePTS <- c(20, 40, 80, 160)
const_vector <- c(rep(0.2, 5), rep(0.7, 15))
Sig <- diag(0.5, d) + matrix(0.5, d, d)
nrep <- 50
# Save SSE results
resu_SSE <- array(0, c(nrep, length(methods), length(NumTimePTS)),
                  dimnames = list(1:nrep, methods, NumTimePTS))
for (idnumT in 1:length(NumTimePTS)) {
  numT = NumTimePTS[idnumT]
  cat("======= T: ", numT, "=======\n")
  for (r in 1:nrep) {
    Ad <- createVARCoefs_ltriangular(p = p, d = d, diag_val = 0.6,
            num_nonzero = d, const_vector = const_vector, range_max = 1)
    Md <- VARparam(Coef = Ad, Sigma = Sig, dof = Inf)
    Y <- simVARmodel(numT = numT, model = Md, burnin = 20)
    EstimRG <- VARshrink(Y, p, 'const', method = 'ridge')
    EstimNS <- VARshrink(Y, p, 'const', method = 'ns')
    EstimFB <- VARshrink(Y, p, 'const', method = 'fbayes', dof = NULL,
                         burnincycle = 1000, mcmccycle = 2000)
    EstimSB <- VARshrink(Y, p, 'const', method = 'sbayes', dof = c(2,5,Inf),
                         prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
    resu_SSE[r, 1, idnumT] <- calcSSE_VARcoef(EstimRG$varparam, Md)
    resu_SSE[r, 2, idnumT] <- calcSSE_VARcoef(EstimNS$varparam, Md)
    resu_SSE[r, 3, idnumT] <- calcSSE_VARcoef(EstimFB$varparam, Md)
    resu_SSE[r, 4, idnumT] <- calcSSE_VARcoef(EstimSB$varparam, Md)
  }
  save(resu_SSE, file = paste0("comparisons_d", d, ".Rdata"))
}
par(mfrow=c(2, 2), mar = c(2.2, 4.2, 2.5, 1))
for (idnumT in 1:length(NumTimePTS)) {
  boxplot(resu_SSE[,, idnumT], log = "y", cex.lab = 1.5, cex.main = 1.5,
        ylab = "SSE", xlab = "", ylim = c(min(resu_SSE), max(resu_SSE)),
        main = paste("T =", NumTimePTS[idnumT]))
  grid()
}

