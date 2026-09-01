# DGP satisfying the control function assumptions: linear first stage and
# jointly normal errors, with Var(U) = 1 marginally.
#   lambda = 0.6, sigma_v = 1  =>  sigma_eps = 0.8
# so the raw (conditional) coefficients are the marginal ones divided by 0.8.
make_cf_data <- function(n = 4000, seed = 5, lambda = 0.6, b = c(0.2, 0.5, -0.4)) {
  set.seed(seed)
  z <- rnorm(n); w <- rnorm(n)
  v <- rnorm(n)
  u <- lambda * v + sqrt(1 - lambda^2) * rnorm(n)   # Var(U) = 1
  x <- 0.9 * z + 0.3 * w + v
  y <- as.numeric(b[1] + b[2] * x + b[3] * w + u >= 0)
  data.frame(y, x, w, z)
}

# Two endogenous regressors. Sigma_v has off-diagonal 0.3 and
# lambda = (0.4, 0.3), so lambda' Sigma_v lambda = 0.322 and
# sigma_eps = sqrt(1 - 0.322).
make_cf_data2 <- function(n = 20000, seed = 21) {
  set.seed(seed)
  Sig <- matrix(c(1, 0.3, 0.3, 1), 2, 2)
  lam <- c(0.4, 0.3)
  q   <- as.numeric(t(lam) %*% Sig %*% lam)
  V   <- matrix(rnorm(2 * n), n, 2) %*% chol(Sig)
  u   <- as.vector(V %*% lam) + sqrt(1 - q) * rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n)
  x1 <- 0.9 * z1 + 0.3 * w + V[, 1]
  x2 <- 0.8 * z2 + 0.2 * w + V[, 2]
  b  <- c(0.2, 0.5, -0.3, -0.4)
  y  <- as.numeric(b[1] + b[2] * x1 + b[3] * x2 + b[4] * w + u >= 0)
  list(data = data.frame(y, x1, x2, w, z1, z2),
       marginal = c(b, lam), sigma_eps = sqrt(1 - q))
}

test_that("the two normalisations differ by exactly the scale factor", {
  d  <- make_cf_data()
  fm <- y ~ x + w | z + w
  cond <- ivcf(fm, d, normalize = "conditional")
  marg <- ivcf(fm, d, normalize = "marginal")

  expect_equal(coef(marg), coef(cond) * marg$scale, tolerance = 1e-12)
  q <- as.numeric(t(marg$lambda) %*% marg$Sigma_v %*% marg$lambda)
  expect_equal(marg$scale, 1 / sqrt(1 + q), tolerance = 1e-12)
  expect_lt(marg$scale, 1)
  expect_equal(cond$scale, 1)
})

test_that("the marginal normalisation recovers the true coefficients", {
  d <- make_cf_data(n = 20000, seed = 12)
  fit <- ivcf(y ~ x + w | z + w, d, normalize = "marginal")
  expect_equal(unname(coef(fit)[c("(Intercept)", "x", "w")]),
               c(0.2, 0.5, -0.4), tolerance = 0.05)
  expect_equal(unname(coef(fit)["V_x"]), 0.6, tolerance = 0.05)
})

test_that("the conditional normalisation recovers beta / sigma_eps", {
  d <- make_cf_data(n = 20000, seed = 12)
  fit <- ivcf(y ~ x + w | z + w, d, normalize = "conditional")
  expect_equal(unname(coef(fit)[c("(Intercept)", "x", "w")]),
               c(0.2, 0.5, -0.4) / 0.8, tolerance = 0.06)
  expect_equal(unname(coef(fit)["V_x"]), 0.75, tolerance = 0.06)
})

test_that("two endogenous regressors are handled", {
  s   <- make_cf_data2()
  fm  <- y ~ x1 + x2 + w | z1 + z2 + w
  marg <- ivcf(fm, s$data, normalize = "marginal")
  cond <- ivcf(fm, s$data, normalize = "conditional")

  expect_equal(marg$k_end, 2L)
  expect_named(coef(marg),
               c("(Intercept)", "x1", "x2", "w", "V_x1", "V_x2"))
  expect_equal(dim(marg$Sigma_v), c(2L, 2L))

  # Marginal normalisation recovers (beta, lambda); conditional recovers them
  # divided by sigma_eps.
  expect_equal(unname(coef(marg)), s$marginal, tolerance = 0.06)
  expect_equal(unname(coef(cond)), s$marginal / s$sigma_eps, tolerance = 0.08)

  # The scale factor is the multivariate one.
  q <- as.numeric(t(cond$lambda) %*% cond$Sigma_v %*% cond$lambda)
  expect_equal(marg$scale, 1 / sqrt(1 + q), tolerance = 1e-12)
  expect_equal(coef(marg), coef(cond) * marg$scale, tolerance = 1e-12)

  expect_true(all(is.finite(diag(vcov(marg)))))
  expect_true(all(diag(vcov(marg)) > 0))
})

test_that("the exogeneity test rejects under endogeneity and not otherwise", {
  s  <- make_cf_data2()
  et <- exogeneity_test(ivcf(y ~ x1 + x2 + w | z1 + z2 + w, s$data))
  expect_equal(et$df, 2L)
  expect_lt(et$p.value, 0.001)

  set.seed(4)
  n <- 5000
  z <- rnorm(n); w <- rnorm(n)
  x <- 0.9 * z + 0.3 * w + rnorm(n)      # exogenous
  y <- as.numeric(0.2 + 0.5 * x - 0.4 * w + rnorm(n) >= 0)
  et0 <- exogeneity_test(ivcf(y ~ x + w | z + w, data.frame(y, x, w, z)))
  expect_equal(et0$df, 1L)
  expect_gt(et0$p.value, 0.01)
})

test_that("with an exogenous regressor the control function reduces to probit", {
  set.seed(9)
  n <- 4000
  z <- rnorm(n); w <- rnorm(n)
  x <- 0.9 * z + 0.3 * w + rnorm(n)
  y <- as.numeric(0.2 + 0.5 * x - 0.4 * w + rnorm(n) >= 0)
  d <- data.frame(y, x, w, z)

  fit <- ivcf(y ~ x + w | z + w, d, normalize = "conditional")
  mle <- glm(y ~ x + w, family = binomial("probit"), data = d)

  expect_lt(abs(fit$lambda), 0.1)
  expect_equal(unname(coef(fit)[c("(Intercept)", "x", "w")]),
               unname(coef(mle)), tolerance = 0.05)
  expect_equal(unname(sqrt(diag(vcov(fit)))[2]),
               unname(sqrt(diag(vcov(mle)))[2]), tolerance = 0.1)
})

test_that("the average structural function is invariant to the normalisation", {
  d  <- make_cf_data()
  fm <- y ~ x + w | z + w
  cond <- ivcf(fm, d, normalize = "conditional")
  marg <- ivcf(fm, d, normalize = "marginal")

  # Under the marginal normalisation Phi(x'beta) IS the ASF.
  expect_equal(predict(marg, type = "response"), predict(marg, type = "asf"),
               tolerance = 1e-12)

  # Under Rivers-Vuong it is not, but rescaling recovers the same ASF.
  expect_equal(predict(cond, type = "asf"), predict(marg, type = "asf"),
               tolerance = 1e-10)
  expect_false(isTRUE(all.equal(predict(cond, type = "response"),
                                predict(cond, type = "asf"),
                                tolerance = 1e-6)))

  # Predictions never involve the first-stage residuals.
  expect_length(predict(marg, type = "response"), nrow(d))
})

test_that("multivariate ASF is also invariant to the normalisation", {
  s  <- make_cf_data2()
  fm <- y ~ x1 + x2 + w | z1 + z2 + w
  expect_equal(predict(ivcf(fm, s$data, normalize = "conditional"), type = "asf"),
               predict(ivcf(fm, s$data, normalize = "marginal"),    type = "asf"),
               tolerance = 1e-10)
})

test_that("marginal effects require the marginal normalisation", {
  skip_if_not_installed("marginaleffects")
  probe <- try(marginaleffects::avg_comparisons(
    glm(am ~ hp, family = binomial, data = mtcars), variables = "hp"),
    silent = TRUE)
  skip_if(inherits(probe, "try-error"),
          "marginaleffects cannot run in this R environment")

  d  <- make_cf_data()
  fm <- y ~ x + w | z + w

  expect_error(
    marginaleffects::avg_slopes(ivcf(fm, d, normalize = "conditional"),
                                variables = "x"),
    "normalize")

  marg <- ivcf(fm, d, normalize = "marginal")
  s <- marginaleffects::avg_slopes(marg, variables = "x")
  expect_true(is.finite(s$estimate))
  expect_true(is.finite(s$std.error))

  # The average structural effect of a continuous regressor is
  # mean(phi(x'beta)) * beta_x.
  b  <- coef(marg)
  bs <- b[names(b) != "V_x"]
  manual <- mean(dnorm(as.vector(marg$parse$X[, names(bs)] %*% bs))) * bs[["x"]]
  expect_equal(s$estimate, manual, tolerance = 1e-5)
})

test_that("ivcf objects behave like other fitted models", {
  d   <- make_cf_data()
  fit <- ivcf(y ~ x + w | z + w, d)
  expect_named(coef(fit), c("(Intercept)", "x", "w", "V_x"))
  expect_equal(names(coef(fit)), rownames(vcov(fit)))
  expect_true(all(diag(vcov(fit)) > 0))
  expect_true(isSymmetric(unname(vcov(fit)), tol = 1e-10))
  expect_equal(nobs(fit), nrow(d))
  expect_output(print(fit), "[Cc]ontrol function")
  expect_output(print(summary(fit)), "Normalisation:")
})
