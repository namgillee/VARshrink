# Example_01.R
# Simulate VAR processes with multivariate t-distribution,
# and run Shrinkage Estimation methods. 
#
# Required packages: MASS, corpcor, ars
#

# Source all files from directories
file.sources = list.files(c("../ShrinkageMethods", "../Functions"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

# Random number generator
set.seed(100)

# Create a VAR model with p=1 and multivariate t-distribution of given dof
myCoef = list(A = list(matrix(c(0.5,0,0,0.5),2,2)), c = c(0.2, 0.7))
myCoef = list(A = list(matrix(c(0.5,0,0,0.5),2,2)), c = c(0.2, 0.7))
myDof = Inf     
myModel = VARparam(Coef = myCoef, Sigma = diag(0.1^2,2), dof = myDof)

# Simulate VAR processes with multivariate t-distribution
lenT = 100
Y = simVARmodel(numT = lenT, model = myModel, burnin = 10)

# Run shrinkage estimation methods
# 1) Multivariate Ridge Regression
EstimRidge = VARshrink(Y, p = 1, const = 'const', method = 'ridge', lambda = NULL)
EstimRidge

# 2) A Nonparametric Shrinkage Method by Opgen-Rein & Strimmer (2007)
EstimNS = VARshrink(Y, p = 1, const = 'const', method = 'ns', lambda = NULL, lambda_var = NULL)
EstimNS

# 3-1) A Full Bayesian Method Using Gibbs MCMC (known fixed dof) by Ni & Sun (2005)
myDof_here = 6
EstimFB1 = VARshrink(Y, p = 1, const = 'const', method = 'fbayes', dof = myDof_here, 
                 burnincycle = 1000, mcmccycle = 2000)
EstimFB1

# 3-2) (unknown dof)
EstimFB2 = VARshrink(Y, p = 1, const = 'const', method = 'fbayes', dof = NULL, 
                 burnincycle = 1000, mcmccycle = 2000)
EstimFB2

# 4-1) A Semiparametric Shrinkage method (known fixed dof) by Lee, Choi & Kim (2016)
myDof_here = 6
EstimSB1 = VARshrink(Y, p = 1, const = 'const', method = 'sbayes', dof = myDof_here, 
                 lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimSB1

# 4-2) (unknown dof)
EstimSB2 = VARshrink(Y, p = 1, const = 'const', method = 'sbayes', dof = NULL, 
                 lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimSB2

# 5) K-fold Cross Validation Alternative for VAR coefficient shrinkage estimation
# 5-1) (known fixed dof)
myDof_here = 6
EstimKCV1 = VARshrink(Y, p = 1, const = 'const', method = 'kcv', dof = myDof_here, 
                  lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimKCV1

# 5-2) (unknown dof)
EstimKCV2 = VARshrink(Y, p = 1, const = 'const', method = 'kcv', dof = NULL, 
                 lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 5, m0 = ncol(Y))
EstimKCV2


print('---- VARshrink: Shrinkage Estimation for VAR models ----')
print('  Results of mean squared errors of estimators     ')
print('')
print(sprintf('Ridge regression       : %.4f', calcMSE_VARcoef(myModel, EstimRidge$varparam)))
print(sprintf('Nonparametric shrinkage: %.4f', calcMSE_VARcoef(myModel, EstimNS$varparam)))
print(sprintf('Full Bayes(fixed dof)  : %.4f', calcMSE_VARcoef(myModel, EstimFB1$varparam)))
print(sprintf('Full Bayes(estim dof)  : %.4f', calcMSE_VARcoef(myModel, EstimFB2$varparam)))
print(sprintf('Semi Bayes(fixed dof)  : %.4f', calcMSE_VARcoef(myModel, EstimSB1$varparam)))
print(sprintf('Semi Bayes(estim dof)  : %.4f', calcMSE_VARcoef(myModel, EstimSB2$varparam)))
print(sprintf('K-fold CV (fixed dof)  : %.4f', calcMSE_VARcoef(myModel, EstimKCV1$varparam)))
print(sprintf('K-fold CV (estim dof)  : %.4f', calcMSE_VARcoef(myModel, EstimKCV2$varparam)))
print('-----------------------------------------------------')



