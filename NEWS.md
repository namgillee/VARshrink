# VARshrink 0.5.0
* Date: 2026-06-17
* [#8](https://github.com/namgillee/VARshrink/pull/8),
* Added conjugate prior to the full Bayesian method, which is available by
  `VARshrink(method = "fbayes", prior_type = "CJ")`.
* Removes linex loss estimates from full Bayesian methods.
* Avoided an unnecessary eigenvalue decomposition computation in full Bayesian
  method: perform the computation only when Q has been updated to not-all-1's.
* Ridge regression estimate has a single lambda value at the minimum of GCV.
* Removed the trivial function `calcSSE_Acoef()`.
* Improved computation of the degrees-of-freedom of shrinkage VAR model.

# VARshrink 0.4.0
* Date: 2026-06-01
* [#7](https://github.com/namgillee/VARshrink/pull/7)
* Move "vars (>= 1.6.1)" from Imports to Depends in `DESCRIPTION`, so that
  inherited methods from **vars** are automatically attached when loading
  **VARshrink**.
* Removed redundant R files due to the update to the **vars** package
  (version 1.6-1): causality_sh.R, restrict_sh.R, roots_sh.R.
* Fixed `VARshrink()` so that the `datamat` element in its output is a
  "data.frame" object. This allows `vars::causality()` to run without errors.
* Added ordinary least squares method to the available `VARshrink()` methods
  with `method = "ols"`.
* Added tools for testing and linting by using **usethis** and **lintr**.
* Updated the full Bayesian estimation method in `lm_full_Bayes_SR()`, so that
  its output additionally includes (1) `$se.param$q`: standard error of
  estimated q, and (2) `$mcmc.param`: MCMC chain of parameters when `store_mcmc=TRUE`.
* Bayesian methods can compute irf confidence bands fast by using an MCMC chain of
  parameters within `irf.varshrinkest()`

# VARshrink 0.3.3
* Date: 2025-11-17
* [#5](https://github.com/namgillee/VARshrink/pull/5)
* [#6](https://github.com/namgillee/VARshrink/pull/6)
* Removed redundant R files due to the update to the **vars** package
  (version 1.6-1): arch.test_sh.R, BQ_sh.R, fevd.varshrinkest.R, h_boot.R,
  h_fecov.R, h_irf.R, normality.test_sh.R, plot.varshirf.R,
  predict.varshrinkest.R.
* Polished the R documentations for all functions.
* README.md includes installation instructions. 
* In the vignette and README.md, `type = "const"` was changed to `type = "none"`
  for `method = "ns"`.
* vignettes/article_varshrink.R was removed.
* LICENSE.txt was removed.

# VARshrink 0.3.2
* Date: 2025-10-13
* [#4](https://github.com/namgillee/VARshrink/pull/4)
* Fixed errors in all methods to allow larger column sizes in the input data
  matrix for season >= 3.
* Fixed `"sbayes"` and `"kcv"` to correctly compute the lag order.
* Fixed `"sbayes"` and `"kcv"` to properly scale the coefficient matrix.

# VARshrink 0.3.1.9100
* Date: 2025-10-07
* [#3](https://github.com/namgillee/VARshrink/pull/3)
* vignettes/article_varshrink.R was added.
* Redundant tokens were removed from .gitignore and .Rbuildignore.
* .Rproj file was removed.

# VARshrink 0.3.1.9000
* Date: 2020-03-09
* Created a release at the GitHub repository for the VARshrink 0.3.1. The version 0.3.1.9000 started.

# VARshrink 0.3.1
* Date: 2019-09-25
* The vignette is reduced in order to compile quicker.
* The R codes include examples.
* In `VARshrink()`, if `method="ns"`, then `type="const"` is switched to
  `type="none"`, and `type="both"` is switched to `type="season"` in order to
  avoid to estimate the constant term. 

# VARshrink 0.3.0
* Date: 2019-08-13
* The pdf article was moved to a separate folder "article_jss/" and ignored from package building.
* The html article was created as a vignette in the "vignettes/" folder.
* Package was built and checked for submission to CRAN.
