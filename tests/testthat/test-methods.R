make_binary_data <- function(n = 1200, seed = 7) {
  set.seed(seed)
  z1 <- rnorm(n); z3 <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x1 <- 0.9 * z1 + 0.4 * z3 + 0.3 * w + 0.6 * u + rnorm(n)
  y  <- as.numeric(0.3 + 0.5 * x1 - 0.4 * w + u >= 0)
  data.frame(y, x1, w, z1, z3)
}

test_that("supported families all fit and return finite standard errors", {
  d <- make_binary_data()
  for (fam in list(binomial("probit"), binomial("logit"), binomial("cloglog"))) {
    fit <- aiv(y ~ x1 + w | z1 + z3 + w, data = d, family = fam)
    expect_true(is.finite(fit$beta_end))
    expect_true(all(is.finite(diag(vcov(fit)))))
    expect_true(all(diag(vcov(fit)) > 0))
  }
})

test_that("coef and vcov are conformable and consistently named", {
  d   <- make_binary_data()
  fit <- aiv(y ~ x1 + w | z1 + z3 + w, data = d)
  expect_equal(names(coef(fit)), rownames(vcov(fit)))
  expect_equal(names(coef(fit)), colnames(vcov(fit)))
  expect_equal(length(coef(fit)), nrow(vcov(fit)))
  expect_true(isSymmetric(unname(vcov(fit)), tol = 1e-10))
})

test_that("predict honours type, newdata and an arbitrary coefficient vector", {
  d   <- make_binary_data()
  fit <- aiv(y ~ x1 + w | z1 + z3 + w, data = d)

  lk <- predict(fit, type = "link")
  rp <- predict(fit, type = "response")
  expect_length(lk, nrow(d))
  expect_true(all(is.finite(lk)))
  expect_true(all(rp > 0 & rp < 1))
  expect_equal(rp, pnorm(lk), tolerance = 1e-12)

  expect_equal(predict(fit, newdata = d[1:10, ], type = "response"),
               rp[1:10], tolerance = 1e-10)

  # Required by the set_coef contract: prediction must respond to an
  # arbitrary coefficient vector, not just the fitted one.
  cf <- coef(fit); cf["x1"] <- cf["x1"] + 1
  expect_equal(predict(fit, type = "link", coefs = cf) - lk, d$x1,
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("print and summary produce glm-style output", {
  d   <- make_binary_data()
  fit <- aiv(y ~ x1 + w | z1 + z3 + w, data = d)

  expect_output(print(fit), "Call:")
  expect_output(print(fit), "Coefficients:")

  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("^Call:", out)))
  expect_true(any(grepl("Estimate\\s+Std\\. Error\\s+z value", out)))
  expect_true(any(grepl("^Endogenous:", out)))
  expect_true(any(grepl("^Identification:", out)))
  # Coefficients must appear exactly once (the earlier version printed the
  # coefficient block twice).
  expect_equal(sum(grepl("^Coefficients:", out)), 1L)
})

test_that("multiple endogenous regressors are supported", {
  set.seed(11)
  n <- 1500
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x1 <- 0.9 * z1 + 0.4 * z3 + 0.3 * w + 0.6 * u + rnorm(n)
  x2 <- 0.8 * z2 + 0.5 * z3 + 0.2 * w + 0.5 * u + rnorm(n)
  y  <- as.numeric(0.3 - 0.4 * w + u >= 0)
  d  <- data.frame(y, x1, x2, w, z1, z2, z3)

  fit <- aiv(y ~ x1 + x2 + w | z1 + z2 + z3 + w, data = d)
  expect_length(fit$beta_end, 2)
  expect_true(all(is.finite(fit$beta_end)))
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_true(fit$overidentified)
})

test_that("the four marginaleffects extension methods honour their contract", {
  # These are the methods we are responsible for. They are tested directly so
  # that contract compliance is verified even where marginaleffects itself
  # cannot run (e.g. a version requiring a newer R than is installed).
  d   <- make_binary_data()
  fit <- aiv(y ~ x1 + w | z1 + z3 + w, data = d)

  expect_equal(auxiv:::get_coef.aiv(fit), coef(fit))
  expect_equal(auxiv:::get_vcov.aiv(fit), vcov(fit))
  expect_null(auxiv:::get_vcov.aiv(fit, vcov = FALSE))

  # set_coef must return an object whose predictions reflect the new values.
  cf  <- coef(fit); cf["x1"] <- cf["x1"] + 1
  fit2 <- auxiv:::set_coef.aiv(fit, cf)
  expect_equal(coef(fit2), cf)
  expect_equal(predict(fit2, type = "link") - predict(fit, type = "link"),
               d$x1, tolerance = 1e-10, ignore_attr = TRUE)

  gp <- auxiv:::get_predict.aiv(fit, newdata = d[1:5, ], type = "response")
  expect_s3_class(gp, "data.frame")
  expect_named(gp, c("rowid", "estimate"))
  expect_equal(nrow(gp), 5L)
  expect_equal(gp$estimate, predict(fit, newdata = d[1:5, ], type = "response"),
               tolerance = 1e-12, ignore_attr = TRUE)

  expect_true("aiv" %in% getOption("marginaleffects_model_classes"))
})

test_that("marginaleffects works and treats a binary treatment as discrete", {
  skip_if_not_installed("marginaleffects")
  # marginaleffects >= 0.20 requires R >= 4.4 (base `%||%`). Skip rather than
  # fail where the installed combination cannot run at all.
  probe <- try({
    dd <- data.frame(a = rbinom(50, 1, 0.5), b = rnorm(50))
    marginaleffects::avg_comparisons(glm(a ~ b, family = binomial, data = dd),
                                     variables = "b")
  }, silent = TRUE)
  skip_if(inherits(probe, "try-error"),
          "marginaleffects cannot run in this R environment")

  set.seed(3)
  n  <- 1500
  z  <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  xs <- 0.9 * z + 0.3 * w + 0.6 * u + rnorm(n)
  tr <- as.numeric(xs >= 0)                    # binary endogenous treatment
  y  <- as.numeric(0.2 + 0.5 * tr - 0.4 * w + u >= 0)
  d  <- data.frame(y, tr, w, z)

  fit <- aiv(y ~ tr + w | z + w, data = d)

  cmp <- marginaleffects::avg_comparisons(fit, variables = "tr")
  expect_true(is.finite(cmp$estimate))
  expect_true(is.finite(cmp$std.error))

  # A 0/1 numeric must be contrasted 0 -> 1, not differentiated. Verify
  # against the discrete difference computed directly from the fit.
  b   <- coef(fit)
  X0  <- X1 <- fit$parse$X
  X0[, "tr"] <- 0; X1[, "tr"] <- 1
  manual <- mean(pnorm(as.vector(X1 %*% b)) - pnorm(as.vector(X0 %*% b)))
  expect_equal(cmp$estimate, manual, tolerance = 1e-6)
})
