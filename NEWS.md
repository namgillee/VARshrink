# VARshrink 0.4.0
* Date: 2025-10-13
* Fixed errors in all methods to allow larger column sizes in the input data
  matrix for season >= 3.
* Fixed errors in sbayes and kcv to correctly compute the lag order and standard
  deviation of the time series, and to properly scale the coefficient matrix.

# VARshrink 0.3.1.9100
* Date: 2025-10-07
* Added vignettes/article_varshrink.R
* Removed redundant tokens from .gitignore and .Rbuildignore
* Removed redundant .Rproj file

# VARshrink 0.3.1.9000
* Date: 2020-03-09
* Created a release at the GitHub repository for the VARshrink 0.3.1. The version 0.3.1.9000 started.

# VARshrink 0.3.1
* Date: 2019-09-25
* The vignette is reduced in order to compile quicker.
* The R codes include examples.
* In VARshrink(), if method="ns", then type="const" is switched to type="none"
and type="both" is switched to type="season" in order to avoid to estimate the constant term. This is because estimated parameters by the 'NS' can result in unexpected errors on the functions borrowed from the package vars.

# VARshrink 0.3.0
* Date: 2019-08-13
* The pdf article was moved to a separate folder "article_jss/" and ignored from package building.
* The html article was created as a vignette in the "vignettes/" folder.
* Package was built and checked for submission to CRAN.
