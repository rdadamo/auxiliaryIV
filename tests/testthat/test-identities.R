# Exact identities. These are the strongest available correctness checks:
# in the linear/Gaussian case the AIV moment reduces to the linear IV moment,
# so the estimator must reproduce IV / 2SLS to machine precision.

make_linear_data <- function(n = 1500, seed = 20260728) {
  set.seed(seed)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x  <- 0.8 * z1 + 0.5 * z2 + 0.4 * w + 0.7 * u + rnorm(n)
  yl <- 1 + 0.5 * x - 0.3 * w + u
  yb <- as.numeric(yl >= 0)
  data.frame(yl, yb, x, w, z1, z2)
}

test_that("just-identified linear AIV reproduces IV exactly", {
  skip_if_not_installed("ivreg")
  d <- make_linear_data()

  fit_iv  <- ivreg::ivreg(yl ~ x + w | z1 + w, data = d)
  fit_aiv <- aiv(yl ~ x + w | z1 + w, data = d, family = gaussian())

  expect_equal(coef(fit_aiv),
               coef(fit_iv)[names(coef(fit_aiv))],
               tolerance = 1e-8)
})

test_that("over-identified linear AIV reproduces 2SLS under the matching weight", {
  skip_if_not_installed("ivreg")
  d <- make_linear_data()
  n <- nrow(d)

  fit_2sls <- ivreg::ivreg(yl ~ x + w | z1 + z2 + w, data = d)

  # The profiled moment is (Z'MZ)^-1 Z'M (y - X_end b), so the effective
  # weight on the 2SLS moment is (Z'MZ)^-1 Omega (Z'MZ)^-1. Choosing
  # Omega = Z'MZ/n collapses it to the 2SLS weight exactly.
  p  <- auxiv:::aiv_parse(yl ~ x + w | z1 + z2 + w, d)
  Mx <- diag(n) - p$X_ex %*% solve(crossprod(p$X_ex), t(p$X_ex))
  Om <- crossprod(p$Z_ex, Mx %*% p$Z_ex) / n

  fit_c <- aiv(yl ~ x + w | z1 + z2 + w, data = d, family = gaussian(),
               Omega = Om)
  expect_equal(coef(fit_c), coef(fit_2sls)[names(coef(fit_c))],
               tolerance = 1e-8)
})

test_that("the default weight does NOT reproduce 2SLS when over-identified", {
  # This is the profiled-vs-full-GMM distinction, and it is a feature, not a
  # bug: the two implementations differ under over-identification.
  skip_if_not_installed("ivreg")
  d <- make_linear_data()
  fit_2sls <- ivreg::ivreg(yl ~ x + w | z1 + z2 + w, data = d)
  fit_c    <- aiv(yl ~ x + w | z1 + z2 + w, data = d, family = gaussian())
  expect_false(isTRUE(all.equal(unname(coef(fit_c)),
                                unname(coef(fit_2sls)[names(coef(fit_c))]),
                                tolerance = 1e-10)))
})

test_that("just-identified probit solves the moment condition", {
  d <- make_linear_data()
  fit <- aiv(yb ~ x + w | z1 + w, data = d, family = binomial("probit"))

  # With one instrument per endogenous regressor the moment equations have a
  # solution, so the sample moment mean(Z_ex * l') vanishes at the optimum.
  p   <- fit$parse
  eta <- as.vector(p$X %*% coef(fit))
  l1  <- auxiv:::index_derivs(fit$family, p$y, eta)$l1
  expect_lt(max(abs(colMeans(p$Z_ex * l1))), 1e-5)

  # Equivalently the objective, which is the weighted norm of the auxiliary
  # coefficients, is driven to zero.
  expect_lt(fit$objective, 1e-10)
})

test_that("the endogenous/exogenous split is inferred correctly", {
  d <- make_linear_data()
  fit <- aiv(yb ~ x + w | z1 + w, data = d, family = binomial("probit"))
  expect_identical(fit$endogenous, "x")
  expect_setequal(fit$exogenous, c("(Intercept)", "w"))
  expect_identical(fit$instruments, "z1")
  expect_false(fit$overidentified)
})

test_that("the order condition is enforced", {
  d <- make_linear_data()
  expect_error(aiv(yb ~ x + w | w, data = d), "order condition")
  expect_error(aiv(yb ~ w | w, data = d), "no endogenous regressors")
})
