# Derivatives of the log-likelihood with respect to the linear index.
#
# Every parameter enters through the index omega = X'beta, so the estimator and
# its asymptotic variance are built entirely from these.

# mu''(eta). Family objects supply mu'(eta) as mu.eta but not this: closed
# forms for the common links, finite differences otherwise.
mu_eta2 <- function(family, eta) {
  switch(family$link,
    probit   = -eta * stats::dnorm(eta),
    logit    = { mu <- family$linkinv(eta); mu * (1 - mu) * (1 - 2 * mu) },
    cloglog  = { e <- exp(eta); exp(eta - e) * (1 - e) },
    log      = exp(eta),
    identity = rep(0, length(eta)),
    {
      h <- pmax(1e-5, 1e-5 * abs(eta))
      (family$mu.eta(eta + h) - family$mu.eta(eta - h)) / (2 * h)
    })
}

# V'(mu). Closed form for the standard families, finite differences otherwise.
variance_deriv <- function(family, mu) {
  switch(family$family,
    binomial = 1 - 2 * mu,
    poisson  = rep(1, length(mu)),
    gaussian = rep(0, length(mu)),
    {
      h <- pmax(1e-6, 1e-6 * abs(mu))
      (family$variance(mu + h) - family$variance(mu - h)) / (2 * h)
    })
}

# First (l1) and second (l2) derivatives of the log-likelihood with respect
# to the index.
index_derivs <- function(family, y, eta) {
  mu  <- family$linkinv(eta)    # mean
  mu1 <- family$mu.eta(eta)     # mu'(eta)
  mu2 <- mu_eta2(family, eta)   # mu''(eta)
  V   <- family$variance(mu)    # variance function: mu(1-mu) binomial, mu Poisson
  r   <- y - mu                 # residual

  list(l1 = r * mu1 / V,
       l2 = -mu1^2 / V + r * (mu2 / V - mu1^2 * variance_deriv(family, mu) / V^2))
}

# Log-likelihood in the index. Used only by the tests, which differentiate it
# numerically to check index_derivs; hence no finite-difference fallback.
index_loglik <- function(family, y, eta) {
  mu <- family$linkinv(eta)
  switch(family$family,
    binomial = y * log(mu) + (1 - y) * log1p(-mu),
    poisson  = y * log(mu) - mu - lgamma(y + 1),
    gaussian = -0.5 * (y - mu)^2,
    stop("unsupported family: ", family$family))
}
