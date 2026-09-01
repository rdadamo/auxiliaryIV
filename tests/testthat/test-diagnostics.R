# Diagnostics that tell the user the estimate is not what it looks like:
# auxiliary coefficients, multiple minima, the edge of the search range, and
# an estimating equation with no root.

healthy <- function(n = 2000, seed = 1) {
  set.seed(seed)
  z <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x <- 0.8 * z + 0.4 * w + 0.7 * u + rnorm(n)
  y <- as.numeric(0.3 - 0.4 * w + u >= 0)
  data.frame(y, x, w, z)
}

# Saturated design: the index is driven far enough that the auxiliary
# coefficients cannot be brought to zero anywhere the model can be fitted.
saturated <- function(n = 5000, b2 = 0.9, seed = 7) {
  set.seed(seed)
  u  <- rnorm(n)
  z  <- (rchisq(n, 10) - 10) / sqrt(20)
  l  <- z + rnorm(n) + u + 2 * (2 * (u > 0) - 1 + u^2 - 1)
  x3 <- rnorm(n) + 0.5 * z^2
  x2 <- l / stats::sd(l)
  x3 <- x3 / stats::sd(x3)
  y  <- as.numeric(1 + b2 * x2 - x3 - u > 0)
  data.frame(y, x2, x3, z)
}

test_that("the auxiliary coefficients are stored and are zero when just identified", {
  f <- aiv(y ~ x + w | z + w, data = healthy())
  expect_named(f$gamma, "z")
  expect_lt(abs(f$gamma[["z"]]), 1e-5)
  expect_lt(f$objective, 1e-10)
})

test_that("a well behaved fit finds one minimum and warns about nothing", {
  expect_silent(f <- aiv(y ~ x + w | z + w, data = healthy()))
  expect_equal(nrow(f$minima), 1L)
  expect_equal(f$minima$beta, unname(f$beta_end))
})

test_that("several local minima are reported and the one closest to zero is returned", {
  # Exercised directly on a bimodal objective: minima at 0.5 and 2.0.
  obj <- function(b) (b - 0.5)^2 * (b - 2)^2
  X   <- matrix(rnorm(100), ncol = 1)
  expect_warning(r <- auxiv:::aiv_outer(obj, X, aiv_control()),
                 "2 local minima")
  expect_equal(nrow(r$minima), 2L)
  expect_equal(sort(round(r$minima$beta, 4)), c(0.5, 2.0))
  expect_equal(r$beta, 0.5, tolerance = 1e-4)
})

# These fits raise several warnings at once — a saturated design also trips
# glm.fit's own. Collect them all rather than asserting on the first.
warnings_from <- function(expr) {
  w <- character()
  withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  w
}

test_that("with distinct minima the global one is reported, without warning", {
  # A shallow trough near 0.5 and the true minimum at 2.0, where the objective
  # is zero. The two are far apart in objective, so the global minimum wins
  # outright and the closest-to-zero tie-break does not apply.
  obj <- function(b) (b - 0.5)^2 * (b - 2)^2 + 0.05 * (b - 2)^2
  X   <- matrix(rnorm(100), ncol = 1)
  w   <- warnings_from(r <- auxiv:::aiv_outer(obj, X, aiv_control()))
  expect_gt(nrow(r$minima), 1L)
  expect_equal(r$beta, 2, tolerance = 1e-4)
  expect_false(any(grepl("same objective", w)))
})

test_that("hitting the edge of the search range warns and names the control", {
  w <- warnings_from(aiv(y ~ x2 + x3 | z + x3, data = saturated(),
                         control = aiv_control(range = 0.3)))
  expect_true(any(grepl("edge of the search range", w)))
  expect_true(any(grepl("aiv_control(range", w, fixed = TRUE)))
})

test_that("stopping at the edge of the fittable region warns, not at the range", {
  # The objective falls without limit but the model stops fitting at 2, well
  # inside the search range. The range is not what binds, so widening it would
  # not help and that warning must not be the one raised.
  obj <- function(b) if (b <= 2) -b else Inf
  X   <- matrix(rnorm(200), ncol = 1)
  w   <- warnings_from(auxiv:::aiv_outer(obj, X, aiv_control()))
  expect_true(any(grepl("edge of the region where the model can be", w)))
  expect_false(any(grepl("edge of the search range", w)))
})

test_that("control accepts a plain list, as glm does", {
  d <- healthy()
  a <- aiv(y ~ x + w | z + w, data = d, control = list(range = 2))
  b <- aiv(y ~ x + w | z + w, data = d, control = aiv_control(range = 2))
  expect_equal(coef(a), coef(b))
  expect_equal(a$algoInfo$range, 2)
})

test_that("no root under just identification warns", {
  w <- warnings_from(aiv(y ~ x2 + x3 | z + x3, data = saturated()))
  expect_true(any(grepl("cannot be driven to zero", w)))
})

test_that("the no-root check also applies with several endogenous regressors", {
  set.seed(4); n <- 4000
  z1 <- rnorm(n); z2 <- rnorm(n); u <- rnorm(n)
  x1 <- 0.8 * z1 + 0.5 * u + rnorm(n)
  x2 <- 0.8 * z2 + 0.5 * u + rnorm(n)
  y  <- as.numeric(0.2 + 0.3 * x1 - 0.3 * x2 + u >= 0)
  d  <- data.frame(y, x1, x2, z1, z2)

  # Two instruments, two endogenous regressors: just identified, so the check
  # is live. It needs the objective at zero, which has no grid to come from in
  # this branch and is evaluated directly.
  f <- aiv(y ~ x1 + x2 | z1 + z2, data = d)
  expect_length(f$beta_end, 2L)
  expect_true(is.finite(f$objective_zero) && f$objective_zero > 0)

  # A root exists here, so it must stay quiet.
  w <- warnings_from(aiv(y ~ x1 + x2 | z1 + z2, data = d))
  expect_false(any(grepl("cannot be driven to zero", w)))
})

test_that("the default optimiser converges tightly enough for the no-root check", {
  # Nelder-Mead stalls once its simplex collapses, leaving the objective around
  # 1e-09 — only a factor of ten below the no-root threshold, close enough to
  # risk a false warning. nlminb reaches about 1e-24 on the same problem, so it
  # is the default, following the gmm package's optfct convention.
  set.seed(4); n <- 4000
  z1 <- rnorm(n); z2 <- rnorm(n); u <- rnorm(n)
  x1 <- 0.8 * z1 + 0.5 * u + rnorm(n)
  x2 <- 0.8 * z2 + 0.5 * u + rnorm(n)
  y  <- as.numeric(0.2 + 0.3 * x1 - 0.3 * x2 + u >= 0)
  d  <- data.frame(y, x1, x2, z1, z2)

  f <- aiv(y ~ x1 + x2 | z1 + z2, data = d)
  expect_identical(f$algoInfo$optimiser, "nlminb")
  expect_identical(f$algoInfo$convergence, 0L)
  expect_lt(f$objective / f$objective_zero, 1e-15)

  # Nelder-Mead remains available and is markedly looser on the same problem.
  g <- aiv(y ~ x1 + x2 | z1 + z2, data = d,
           control = aiv_control(optfct = "optim"))
  expect_identical(g$algoInfo$optimiser, "optim (Nelder-Mead)")
  expect_gt(g$objective / g$objective_zero, f$objective / f$objective_zero)
  expect_equal(unname(coef(g)), unname(coef(f)), tolerance = 1e-3)
})

test_that("the grid branch records how the search went", {
  f <- aiv(y ~ x + w | z + w, data = healthy())
  a <- f$algoInfo
  expect_identical(a$optimiser, "grid search + optimize")
  expect_identical(a$range, 4)
  expect_gt(a$admissible, 0L)
  expect_lte(a$admissible, a$counts)
  expect_identical(a$local_minima, 1L)

  # The half-width is the range in standard deviations of the endogenous
  # regressor, and the region where the model actually fits sits inside it.
  expect_equal(unname(a$half), 4 / sd(healthy()$x), tolerance = 1e-8)
  expect_gte(a$admissible_range[1], -a$half)
  expect_lte(a$admissible_range[2],  a$half)
})

test_that("a positive objective under over-identification does not warn", {
  # Two instruments for one endogenous regressor: gamma = 0 is two equations
  # in one unknown, so a strictly positive objective is expected, not a fault.
  set.seed(3); n <- 2000
  z1 <- rnorm(n); z2 <- rnorm(n); w <- rnorm(n); u <- rnorm(n)
  x <- 0.6 * z1 + 0.5 * z2 + 0.4 * w + 0.7 * u + rnorm(n)
  y <- as.numeric(0.3 - 0.4 * w + u >= 0)
  d <- data.frame(y, x, w, z1, z2)
  expect_silent(f <- aiv(y ~ x + w | z1 + z2 + w, data = d))
  expect_true(f$overidentified)
  expect_gt(f$objective, 0)
})
