# Source all files from directories
file.sources = list.files(c("../ShrinkageMethods", "../Functions"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

# Random number generator
set.seed(100)

p <- 1
d <- 40
numT <- 80
nrep <- 50
methods <- c("ridge", "ns", "fbayes", "sbayes")
const_vector<-c(rep(0.2,15),rep(0.8,25))
Sig <- diag(0.5,40)+matrix(0.5,40,40)
dof <- Inf



# 결과 저장 
resu_SSE <- matrix(0, nrow = nrep, ncol = length(methods))


for (r in 1:nrep) {
  print(paste("Repeat:", r))

  A40 <- createVARCoefs_ltriangular(p = p,d = d,num_nonzero = d,const_vector = const_vector )
  M40 <- VARparam(Coef = A40, Sigma = Sig, dof = dof)
  V40 <- simVARmodel(numT = numT,model = M40, burnin = 20)
  V40_ridge <- VARshrink(Y = V40,p = p,const_type = 'const',method = 'ridge',lambda = NULL)
  V40_ns <- VARshrink(Y = V40,p = p,const_type = 'const',method = 'ns',lambda = NULL,lambda_var = NULL)
  V40_fbayes <- VARshrink(Y = V40,p = p,const_type = 'const',method = 'fbayes',dof = NULL, 
                          burnincycle = 1000, mcmccycle = 2000)
  V40_sbayes <- VARshrink(Y = V40,p = p,const_type = 'const',method = 'sbayes', dof = NULL, 
                          lambda = NULL, lambda_var = NULL, prior_type = 'NCJ', num_folds = 5, m0 = ncol(V40))
  resu_SSE[r,1] <- calcSSE_VARcoef(model1 = V40_ridge$varparam, model2 = M40, include_const_vector = TRUE)
  resu_SSE[r,2] <- calcSSE_VARcoef(model1 = V40_ns$varparam, model2 = M40, include_const_vector = TRUE)
  resu_SSE[r,3] <- calcSSE_VARcoef(model1 = V40_fbayes$varparam, model2 = M40, include_const_vector = TRUE)
  resu_SSE[r,4] <- calcSSE_VARcoef(model1 = V40_sbayes$varparam, model2 = M40, include_const_vector = TRUE)
}

save(resu_SSE, file = "comparisons_d40T80.Rdata")

colnames(resu_SSE)<-methods
boxplot(resu_SSE,ylab="MSE")