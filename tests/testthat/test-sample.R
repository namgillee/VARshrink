test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Bayesian IRF returns vars-compatible varirf object", {
  data(Canada, package = "vars")
  fit <- VARshrink(diff(Canada), p = 1, method = "ridge")
  coef_draw <- t(Bcoef_sh(fit))
  expect_gt(nrow(coef_draw), fit$K * fit$p)
  sigma_draw <- crossprod(resid(fit)) / df.residual(fit$varresult[[1]])
  fit$mcmc.param <- list(
    list(Psi = coef_draw, Sigma = sigma_draw),
    list(Psi = coef_draw * 1.01, Sigma = sigma_draw * 1.02),
    list(Psi = coef_draw * 0.99, Sigma = sigma_draw * 0.98)
  )

  result <- irf(fit, impulse = "e", response = c("prod", "rw"),
                n.ahead = 3, boot = TRUE, ci = 0.9)

  expect_s3_class(result, "varirf")
  expect_named(result$irf, "e")
  expect_equal(dim(result$irf$e), c(4, 2))
  expect_equal(dim(result$Lower$e), c(4, 2))
  expect_equal(dim(result$Upper$e), c(4, 2))
  expect_equal(result$runs, 3)
  expect_equal(result$ci, 0.1)
})

test_that("Bayesian IRF supports non-orthogonal point IRF without bands", {
  data(Canada, package = "vars")
  fit <- VARshrink(diff(Canada), p = 1, method = "ridge")
  fit$mcmc.param <- list(list(Psi = t(Bcoef_sh(fit))))

  result <- irf(fit, impulse = "e", response = "prod", n.ahead = 0,
                ortho = FALSE, boot = FALSE)

  expect_s3_class(result, "varirf")
  expect_equal(dim(result$irf$e), c(1, 1))
  expect_null(result$Lower)
  expect_null(result$Upper)
})
