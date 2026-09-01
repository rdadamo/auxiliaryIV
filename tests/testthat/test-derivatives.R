test_that("index derivatives match numerical differentiation", {
  skip_if_not_installed("numDeriv")

  fams <- list(binomial("probit"), binomial("logit"), binomial("cloglog"),
               poisson("log"), gaussian("identity"))

  for (fam in fams) {
    eta <- seq(-2.5, 2.5, length.out = 21)
    yv <- if (fam$family == "binomial") rep(c(0, 1), length.out = length(eta))
          else if (fam$family == "poisson") rep(c(0, 1, 2, 5), length.out = length(eta))
          else seq(-1, 1, length.out = length(eta))

    d <- auxiv:::index_derivs(fam, yv, eta)

    num1 <- vapply(seq_along(eta), function(i)
      numDeriv::grad(function(e) auxiv:::index_loglik(fam, yv[i], e), eta[i]),
      numeric(1))
    num2 <- vapply(seq_along(eta), function(i)
      numDeriv::hessian(function(e) auxiv:::index_loglik(fam, yv[i], e),
                        eta[i])[1, 1], numeric(1))

    expect_equal(d$l1, num1, tolerance = 1e-5,
                 info = paste("l' for", fam$family, fam$link))
    expect_equal(d$l2, num2, tolerance = 1e-4,
                 info = paste("l'' for", fam$family, fam$link))
  }
})

test_that("probit derivatives match the closed form used in the paper", {
  fam <- binomial("probit")
  eta <- seq(-3, 3, length.out = 41)
  yv  <- rep(c(0, 1), length.out = length(eta))

  F_ <- pnorm(eta); F1 <- dnorm(eta); F2 <- -eta * dnorm(eta)
  legacy_l1 <- yv * F1 / F_ - (1 - yv) * F1 / (1 - F_)
  legacy_l2 <- yv * (F2 / F_ - F1^2 / F_^2) -
               (1 - yv) * (F2 / (1 - F_) + F1^2 / (1 - F_)^2)

  d <- auxiv:::index_derivs(fam, yv, eta)
  expect_equal(d$l1, legacy_l1, tolerance = 1e-12)
  expect_equal(d$l2, legacy_l2, tolerance = 1e-12)
})

test_that("analytic mu''(eta) matches numerical differentiation", {
  for (fam in list(binomial("probit"), binomial("logit"),
                   binomial("cloglog"), poisson("log"))) {
    eta <- seq(-2, 2, length.out = 21)
    h <- pmax(1e-5, 1e-5 * abs(eta))
    num <- (fam$mu.eta(eta + h) - fam$mu.eta(eta - h)) / (2 * h)
    expect_equal(auxiv:::mu_eta2(fam, eta), num, tolerance = 1e-4,
                 info = fam$link)
  }
})
