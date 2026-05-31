test_that("Bayesian IRF returns vars-compatible varirf object", {
  data(Canada, package = "vars")
  fit <- VARshrink(diff(Canada), p = 1, method = "fbayes", type = "const",
                   burnincycle = 20, mcmccycle = 50, store_mcmc = TRUE)

  # Bayesian irf with CI
  result_bayes <- irf(fit, impulse = "e", response = c("prod", "rw"),
                n.ahead = 3, ortho = FALSE, boot = TRUE, ci = 0.9)

  expect_s3_class(result_bayes, "varirf")
  expect_named(result_bayes$irf, "e")
  expect_equal(dim(result_bayes$irf$e), c(4, 2))
  expect_equal(dim(result_bayes$Lower$e), c(4, 2))
  expect_equal(dim(result_bayes$Upper$e), c(4, 2))
  expect_equal(result_bayes$runs, 50)
  expect_equal(result_bayes$ci, 0.1)

  # Compare point irfs of Bayesian irf and vars::irf when ortho=FALSE
  fit_nomcmc <- fit
  fit_nomcmc$mcmc.param <- NULL
  result_nomcmc <- irf(fit_nomcmc, impulse = "e", response = c("prod", "rw"),
                     n.ahead = 3, ortho = FALSE, boot = FALSE)

  expect_equal(result_bayes$irf$e, result_nomcmc$irf$e)
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
