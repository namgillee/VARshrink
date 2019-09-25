#' Create coefficients of a VAR model
#'
#' Randomly create sparse lower-triangular matrices for VAR coefficients of
#' lagged endogenous variables, and set a constant vector.
#'
#' Consider VAR(p) model:
#' \deqn{y_t = A_1 y_{t-1} + ... + A_p y_{t-p} + c + e_t,}
#' with the constant deterministic variable (d_t = 1).
#' The function creates the coefficient matrices A_1, ..., A_p and constant
#' vector c.
#'
#' Diagonal elements of each K-by-K matrix A_k are all equal to diag_val,
#' and off-diagonal elements are all zero except for a few randomly selected
#' nonzero elements. Nonzero off-diagonal elements are selected from
#' lower-triangular parts of A_i and the values are drawn from a uniform
#' distribution over [-range_max, -range_min] U [range_min, range_max].
#'
#' @param p lag order
#' @param K Number of time series variables.
#' @param diag_val diagonal values of A1,...,Ap
#' @param num_nonzero Number of nonzero entries on the lower-triangular parts of
#' A1, ..., Ap
#' @param const_vector constant vector c of the VAR model
#' @param range_min,range_max Each nonzero off-diagonal entry of coefficient
#' matrices is drawn uniformly from the interval
#' [-range_max, -range_min] U [range_min, range_max]
#' @return A list object with components $A and $c. $A is a list of K-by-K
#' matrices A_1, ..., A_p, and $c is a constant vector of length K.
#' @importFrom stats runif
#' @examples
#' p <- 1; K <- 20;
#' const_vector <- c(rep(0.2, 5), rep(0.7, 15))
#' createVARCoefs_ltriangular(p = p, K = K, diag_val = 0.6,
#' num_nonzero = K, const_vector = const_vector, range_max = 1)
#' @export
createVARCoefs_ltriangular <- function(p = 1, K = 5, diag_val = 1 / p,
                                       num_nonzero = 0, const_vector = NULL,
                                       range_min = 0.2, range_max = 1 / p) {
  #-- Generate coefficients ---#
  var_coef <- vector("list", 2)
  names(var_coef) <- c("A", "c")

  coef <- matrix(0, K, K * p)
	for (j in 1:p) {
    diag(coef[, (1 + (j - 1) * K):(j * K)]) <- diag_val
  }

  #-- Draw nonzero coefficients ---%
  if (floor(num_nonzero) > 0) {
    #INDEX FOR LOWER-TRIANGULAR PARTS OF EACH A_i,i=1,..,p
    idxzero <- NULL
    lower_indices <- which(lower.tri(matrix(0, K, K)))
    for (j in 1:p) {
      idxzero <- c(idxzero, lower_indices + (j - 1) * K * K)
    }
    lenzero <- length(idxzero)

    idx_nonzero <- sample(lenzero,  min(lenzero, floor(num_nonzero)))
  	val_nonzero <- runif(length(idx_nonzero), min = -1, max = 1)
  	val_nonzero[val_nonzero >= 0] <- (range_max - range_min) *
  	  val_nonzero[val_nonzero >= 0] + range_min
  	val_nonzero[val_nonzero < 0] <- (range_max - range_min) *
  	  val_nonzero[val_nonzero < 0] - range_min

  	coef[idxzero[idx_nonzero]] <- val_nonzero
  }

  #-- Return value --#
	for (j in 1:p) {
		var_coef$A[[j]] <- coef[, (1 + (j - 1) * K):(j * K)]
	}

	if (is.null(const_vector)) {
		var_coef$c <- matrix(0, K, 1)
	}	else 	{
		var_coef$c <- const_vector
	}

	return(var_coef)
}
