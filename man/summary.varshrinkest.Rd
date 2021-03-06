% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.varshrinkest.R
\name{summary.varshrinkest}
\alias{summary.varshrinkest}
\title{Summary method for an object of class 'varshrinkest',
VAR parameters estimated by VARshrink()}
\usage{
\method{summary}{varshrinkest}(object, equations = NULL, ...)
}
\arguments{
\item{object}{An object of class "varshrinkest", usually
a result of call to "VARshrink()".}

\item{equations}{Subset of names of endogenous time series variables
to summarize.}

\item{...}{Currently not used.}
}
\description{
Extend summary.varest() to class 'varshrinest' to incorporate
adapted methods for new classes:
summary.shrinklm(), logLik.varshrinkest(), roots.varshrinkest().
}
\details{
Code is modified to avoid call to data matrices ($y, $datamat)
and to use effective numbers of parameters of shrinkage estimates.

Output includes the scale matrix, Sigma, and degree-of-freedom, dof,
for multivariate t-distribution for residuals.
}
\examples{
data(Canada, package = "vars")
y <- diff(Canada)
estim <- VARshrink(y, p = 2, type = "const", method = "ridge")
summary(estim)
}
