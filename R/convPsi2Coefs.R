#' @export
convPsi2Coefs <- function(Psi)
{
  # Convert (or break down) a coefficient matrix Psi
  #into a list of $A and $c.
  #In specific,
  #         t(Psi) = cbind(c,A[[1]],...,A[[p]])
  # Input
  #    Psi    : a (1+dp)xd matrix
  # Output
  #    Coef   : $A is a list of length p,  and $c is a vector of length d
  #

  s = dim(Psi)
  d = s[2]
  p = s[1] %/% d

  Coefs = vector('list',2)
  names(Coefs) = c('A','c')

  is_const = s[1] %% d  #if ==1, there is constant vector c.
  if (is_const == 0) {  #if ==0, there is no constant vector c.
    f = rep( 1:p, rep(d^2,p) )
    PsiSplit = split(t(Psi), f = f)
    Coefs$A = lapply(PsiSplit, matrix, d)
    Coefs$c = NULL
  } else {
    f = rep( 1:(p+1), c(d, rep(d^2,p)) )
    PsiSplit = split(t(Psi), f = f)
    Coefs$A = lapply(PsiSplit[-1], matrix, d)
    Coefs$c = matrix(PsiSplit[[1]], d)
  }

  return(Coefs)
}
