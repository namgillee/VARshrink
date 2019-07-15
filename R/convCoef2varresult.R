# coefficients	a named vector of coefficients
# residuals	the residuals, that is response minus fitted values.
# fitted.values	the fitted mean values.
# rank		the numeric rank of the fitted linear model.
# df.residual	the residual degrees of freedom.
# call		the matched call.
#



convCoef2varresult <- function(Coef,Y,df.residual,call=NULL)
{
  d <- ncol(Y)
  p <- length(Coef$A)

  # Psi0 = cbind(A1, A2, ..., Ap, c), or, cbind(A1, A2, ..., Ap)
  Psi0 <- ifelse(is.null(Coef$c), matrix(0,d,d*p), matrix(0,d,d*p+1))
  for(j in 1:p){

    Psi0[,(d*(j-1)+1):(d*j)] <- Coef$A[[j]]

    }
  if (!is.null(Coef$c)) {
    Psi0[,d*p+1] <- Coef$c
  }

  Num <- vector(0,d*p)
  for(z in 1:p){

    Num[(d*(z-1)+1):(d*z)] <- rep(z,d)

  }
  colnames(Psi0) <- paste0(colnames(Y),".l",Num)

  #datY, datX
  datY = Y[-(1:p),]
  datX  = NULL
  for (h in 1:p) {
    datX = cbind(datX, Y[(p+1-h):(nrow(Y)-h),])
  }
  if (!is.null(Coef$c)) {
    # Append 1 to datX, because we will estimate
    # constant vector by shrinkage
    datX = cbind(datX, 1)
  }

  residuals <- datY-datX%*%Psi0

  fitted <- datX%*%Psi0


  varresult = vector('list', d)
  names(varresult) = ifelse(!is.null(colnames(Y)), colnames(Y), paste("Y",1:d,sep = ""))



  for(i in 1:d){

    varresult[[i]]$coefficients <- Psi0[i,]
    varresult[[i]]$residuals <- resuduals[i,]
    varresult[[i]]$fitted.values <- fitted[i,]
    varresult[[i]]$rank <- ncol(Y)*p
    varresult[[i]]$df.residual <- df.residual
    varresult[[i]]$call <- call #list

    class(varresult[[i]]) <- c('lm')

  }



  return(varresult)
}
