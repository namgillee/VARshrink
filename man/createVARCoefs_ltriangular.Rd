% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createVARCoefs_ltriangular.R
\name{createVARCoefs_ltriangular}
\alias{createVARCoefs_ltriangular}
\title{Create coefficients of a VAR model}
\usage{
createVARCoefs_ltriangular(p = 1, K = 5, diag_val = 1/p,
  num_nonzero = 0, const_vector = NULL, range_min = 0.2,
  range_max = 1/p)
}
\arguments{
\item{p}{lag order}

\item{K}{Number of time series variables.}

\item{diag_val}{diagonal values of A1,...,Ap}

\item{num_nonzero}{Number of nonzero entries on the lower-triangular parts of
A1, ..., Ap}

\item{const_vector}{constant vector c of the VAR model}

\item{range_min, range_max}{Each nonzero off-diagonal entry of coefficient
matrices is drawn uniformly from the interval
[-range_max, -range_min] U [range_min, range_max]}
}
\value{
A list object with components $A and $c. $A is a list of K-by-K
matrices A_1, ..., A_p, and $c is a constant vector of length K.
}
\description{
Randomly create sparse lower-triangular matrices for VAR coefficients of
lagged endogenous variables, and set a constant vector.
}
\details{
Consider VAR(p) model:
\deqn{y_t = A_1 y_{t-1} + ... + A_p y_{t-p} + c + e_t,}
with the constant deterministic variable (d_t = 1).
The function creates the coefficient matrices A_1, ..., A_p and constant
vector c.

Diagonal elements of each K-by-K matrix A_k are all equal to diag_val,
and off-diagonal elements are all zero except for a few randomly selected
nonzero elements. Nonzero off-diagonal elements are selected from
lower-triangular parts of A_i and the values are drawn from a uniform
distribution over [-range_max, -range_min] U [range_min, range_max].
}
\examples{
p <- 1; K <- 20;
const_vector <- c(rep(0.2, 5), rep(0.7, 15))
createVARCoefs_ltriangular(p = p, K = K, diag_val = 0.6,
num_nonzero = K, const_vector = const_vector, range_max = 1)
}
