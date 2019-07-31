#' A variant of vars::Bcoef() for an object of class 'varshrinkest'
#'
#' Code is modified to avoid call to data matrices ($y, $datamat)
#' and to use effective numbers of parameters of shrinkage estimates.

Bcoef.varshrinkest <- function (x) {

  if (!inherits(x, "varest")) {
    stop("\nPlease provide an object inheriting class 'varest'.\n")
  }
  y.names <- names(x$varresult)
  x.names <- if (x$K > 0) {
    names(coef(x$varresult[[1]]))
  } else {
    character(0)
  }
  B <- matrix(0, nrow = x$K, ncol = length(x.names))
  if (is.null(x$restriction)) {
    for (i in 1:x$K) {
      B[i, ] <- coef(x$varresult[[i]])
    }
  }
  else if (!(is.null(x$restriction))) {
    for (i in 1:x$K) {
      restrictions <- x$restrictions
      restrictions[i, restrictions[i, ] == TRUE] <- coef(x$varresult[[i]])
      temp <- restrictions[i, ]
      B[i, ] <- temp
    }
  }
  colnames(B) <- x.names
  rownames(B) <- y.names
  return(B)
}
