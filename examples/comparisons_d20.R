# Source all files from directories
file.sources = list.files(c("../ShrinkageMethods", "../Functions"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
sapply(file.sources,source,.GlobalEnv)

# Random number generator
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


