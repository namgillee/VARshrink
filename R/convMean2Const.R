#' @export
convMean2Const <- function(ybar, myCoefs)
{
  # Convert the mean vector ('ybar') into
  #a constant vector ('const') for VAR model
  #         const = (I - sum_i(A_i)) * ybar

  d = ncol(myCoefs$A[[1]])
  p = length(myCoefs$A)

  A = diag(d)
  for (i in 1:p)
    A = A - myCoefs$A[[i]]

  const = A %*% ybar
  return(const)
}
