# VARshrink: Shrinkage Estimation Methods for Vector Autoregressive (VAR) Models. 

VAR model is a fundamental and effective approach for multivariate time series analysis.
Shrinkage estimation methods can be applied to high-dimensional VAR models with dimensionality greater than the number of observations, contrary to the standard ordinary least squares method. 

The package **VARshrink** aims to be an integrative R package delivering nonparametric, parametric, and semiparametric methods in a unified and consistent manner.

The package **VARshrink** provides a simple interface function `VARshrink()`, which is an extension of the function `VAR()` in the **vars** package. 

Example:

```
data(Canada, package = "vars")
y <- diff(Canada)
estim <- VARshrink(Y, p = 2, type = "const", method = "ns")
plot(predict(estim), names = "U")
```

### References

N. Lee, H. Choi, and S.-H. Kim (2016). Bayes shrinkage estimation for high-dimensional VAR models with scale mixture of normal distributions for noise. Computational Statistics & Data Analysis 101, 250-276. doi: 10.1016/j.csda.2016.03.007

S. Ni and D. Sun (2005). Bayesian estimates for vector autoregressive models. Journal of Business & Economic Statistics 23(1), 105-117. doi: 10.1198/073500104000000622

R. Opgen-Rhein and K. Strimmer (2007). Learning causal networks from systems biology time
course data: an effective model selection procedure for the vector autoregressive process.
BMC Bioinformatics 8(2), S3. doi: 10.1186/1471-2105-8-S2-S3.
