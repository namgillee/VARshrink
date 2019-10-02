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
