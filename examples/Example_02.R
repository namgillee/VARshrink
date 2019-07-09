# Source all files from directories
file.sources = list.files(c("../ShrinkageMethods", "../Functions"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

set.seed(100)
#Canada
data(Canada, package = "vars")

#const = 'const', method = 'ridge'
EstimRidge2 = VARshrink(Canada, p = 2, const = 'const', method = 'ridge')
EstimRidge2

#const = 'const', method = 'ns'
Estimns2 = VARshrink(Canada, p = 2, const = 'const', method = 'ns')
Estimns2

EstimFB2 = VARshrink(Canada, p = 2, const = 'const', method = 'fbayes', dof = NULL, 
                     burnincycle = 1000, mcmccycle = 2000)
EstimFB2

EstimSB2 = VARshrink(Canada, p = 2, const = 'const', method = 'sbayes', dof = NULL, 
                     lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 10, m0 = ncol(Canada))
EstimSB2


# #BVAR
# data("mts-examples",package="MTS")
# z=log(qgdp[,3:5])
# zt=diffM(z)*100
# C=0.1*diag(rep(1,7))
# V0=diag(rep(1,3))
# BVAR(zt,p=2,C,V0)

# myDof_here = 6
# EstimFB1 = VARshrink(zt, p = 1, const = 'const', method = 'fbayes', dof = myDof_here, 
#                      burnincycle = 1000, mcmccycle = 2000)
# EstimFB1
# 
# EstimFB2 = VARshrink(zt, p = 1, const = 'const', method = 'fbayes', dof = NULL, 
#                      burnincycle = 1000, mcmccycle = 2000)
# EstimFB2
# 
# EstimSB2 = VARshrink(zt, p = 1, const = 'const', method = 'sbayes', dof = NULL, 
#                      lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 10, m0 = ncol(Canada))
# EstimSB2


