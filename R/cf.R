# Control function estimator (Rivers and Vuong, 1988), for comparison with
# aiv().
#
#   y = 1{X'beta + U >= 0},   x_end = Pi'Z + V,   (U, V) jointly normal.
#
# Writing U = lambda'V + eps with eps independent of V,
#   P(y = 1 | X, V) = Phi((X'beta + lambda'V) / sigma_eps),
# so a probit of y on (X, Vhat) estimates beta/sigma_eps and lambda/sigma_eps.
# Which variance is normalised to one therefore fixes the scale:
#
#   "conditional"  Var(U | Z, V) = 1; the probit coefficients are the estimates
#   "marginal"     Var(U) = 1; every coefficient is multiplied by
#                  sigma_eps = (1 + lambda'Sigma_v lambda)^(-1/2)
#
# Standard errors correct for estimation of the first stage, and under the
# marginal normalisation also for the estimation of Sigma_v.

#' Control function estimator for a binary outcome
#'
#' Two-step control function (Rivers and Vuong, 1988): the first-stage
#' residuals are added to a probit as extra regressors. Any number of
#' continuous endogenous regressors is supported.
#'
#' @param formula A two-part formula, `y ~ x_end + x_ex | z_ex + x_ex`, with
#'   the complete instrument set in the second part, as in [aiv()].
#' @param data A data frame.
#' @param normalize Which variance of the structural error is set to one.
#'   `"marginal"` (the default) sets \eqn{Var(U) = 1}, matching [aiv()];
#'   `"conditional"` sets \eqn{Var(U \mid Z, V) = 1}, as in Rivers and Vuong.
#' @return An object of class `"ivcf"`. The coefficients on the first-stage
#'   residuals are named `V_<regressor>`; a test of exogeneity is a test that
#'   they are jointly zero.
#' @references Rivers, D. and Q. H. Vuong (1988). Limited information
#'   estimators and exogeneity tests for simultaneous probit models.
#'   \emph{Journal of Econometrics} 39, 347--366.
#' @examples
#' set.seed(1)
#' n <- 1000
#' z <- rnorm(n); w <- rnorm(n)
#' v <- rnorm(n); u <- 0.6 * v + sqrt(1 - 0.36) * rnorm(n)
#' x <- 0.9 * z + 0.3 * w + v
#' y <- as.numeric(0.2 + 0.5 * x - 0.4 * w + u >= 0)
#' summary(ivcf(y ~ x + w | z + w, data = data.frame(y, x, w, z)))
#' @export
ivcf <- function(formula, data, normalize = c("marginal", "conditional")) {
  normalize <- match.arg(normalize)

  p  <- aiv_parse(formula, data)
  ke <- ncol(p$X_end)
  y  <- p$y
  n  <- length(y)
  X  <- p$X                      # structural regressors
  Zf <- cbind(p$Z_ex, p$X_ex)    # complete instrument set

  # First stage: one regression per endogenous regressor. With identical
  # regressors across equations, equation-by-equation least squares is
  # efficient.
  Pi <- qr.solve(Zf, p$X_end)
  V  <- p$X_end - Zf %*% Pi
  colnames(V) <- paste0("V_", colnames(p$X_end))

  # Second stage: probit of y on (X, V).
  W   <- cbind(X, V)
  fam <- stats::binomial("probit")
  fit <- stats::glm.fit(x = W, y = y, family = fam, intercept = FALSE)
  theta <- fit$coefficients
  names(theta) <- colnames(W)

  d      <- index_derivs(fam, y, as.vector(W %*% theta))
  is_v   <- colnames(W) %in% colnames(V)
  lambda <- theta[is_v]

  # Second-stage Hessian and the first-stage correction. Differentiating the
  # second-stage score with respect to vec(Pi) gives lambda' %x% M1 with
  # M1 = -E[l'' W Z'], so the per-observation correction reduces to the scalar
  # lambda'V_i times M1 Q^{-1} Z_i. It is proportional to lambda and so
  # vanishes when the regressors are exogenous.
  M_theta <- crossprod(W * d$l2, W) / n
  M1      <- -crossprod(W * d$l2, Zf) / n
  Q       <- crossprod(Zf) / n
  Vlam    <- as.vector(V %*% lambda)

  u <- W * d$l1 + (Zf * Vlam) %*% solve(Q, t(M1))
  M_theta_inv <- solve(M_theta)
  V_theta <- M_theta_inv %*% (crossprod(u) / n) %*% M_theta_inv / n
  dimnames(V_theta) <- list(names(theta), names(theta))

  Sigma_v <- crossprod(V) / n
  q       <- as.numeric(t(lambda) %*% Sigma_v %*% lambda)
  sc      <- 1 / sqrt(1 + q)

  if (normalize == "conditional") {
    est <- theta
    vc  <- V_theta
  } else {
    # Joint variance of (theta, vec(Sigma_v)). The influence function of
    # Sigma_v is vec(V_i V_i' - Sigma_v); the first-stage correction to it is
    # of smaller order and drops out.
    j_idx <- rep(seq_len(ke), times = ke)
    l_idx <- rep(seq_len(ke), each  = ke)
    psi_S <- V[, j_idx, drop = FALSE] * V[, l_idx, drop = FALSE] -
      matrix(as.vector(Sigma_v), n, ke^2, byrow = TRUE)

    psi  <- cbind(-u %*% M_theta_inv, psi_S)
    V_xi <- crossprod(psi) / n^2

    # d(sc)/d(theta) is nonzero only in the lambda positions, and
    # d(sc)/d(vec(Sigma_v)) = -sc^3 vec(lambda lambda') / 2.
    dc_dtheta       <- numeric(length(theta))
    dc_dtheta[is_v] <- -sc^3 * as.vector(Sigma_v %*% lambda)
    dc_dvecS        <- -sc^3 / 2 * as.vector(tcrossprod(lambda))

    J <- cbind(sc * diag(length(theta)) + outer(theta, dc_dtheta),
               outer(theta, dc_dvecS))

    est <- theta * sc
    vc  <- J %*% V_xi %*% t(J)
    dimnames(vc) <- list(names(theta), names(theta))
  }

  structure(list(
    coefficients = est, vcov = vc,
    theta = theta, vcov_conditional = V_theta,
    lambda = lambda, Sigma_v = Sigma_v,
    scale = if (normalize == "marginal") sc else 1,
    normalize = normalize, residuals_first_stage = V, Pi = Pi,
    family = fam, formula = formula, terms = p$terms, call = match.call(),
    endogenous = p$end, exogenous = p$ex, instruments = p$iv,
    n = n, k_end = ke, parse = p
  ), class = "ivcf")
}

#' @export
coef.ivcf <- function(object, ...) object$coefficients

#' @export
vcov.ivcf <- function(object, ...) object$vcov

#' @export
nobs.ivcf <- function(object, ...) object$n

#' Predictions from a control function fit
#'
#' The quantity of interest for marginal effects is the average structural
#' function \eqn{ASF(x) = E_U[1\{x'\beta + U \ge 0\}] = \Phi(x'\beta)}, which
#' requires \eqn{Var(U) = 1} \emph{marginally}. Under the marginal
#' normalisation \eqn{\Phi(x'\beta)} is therefore the ASF and `"response"` and
#' `"asf"` coincide.
#'
#' Under the Rivers-Vuong normalisation they do not. There
#' \eqn{\Phi(x'\beta)} is \eqn{P(y = 1 \mid X, V = 0)}, not a structural
#' quantity. Integrating out the first-stage error recovers the ASF,
#' \deqn{E_V[\Phi(x'\beta + \lambda'V)] = \Phi(x'\beta / \sqrt{1 +
#'       \lambda'\Sigma_v\lambda}),}
#' so `"asf"` applies that rescaling. The ASF itself is invariant to the
#' normalisation; only the scale of the coefficients differs.
#'
#' Predictions never involve the first-stage residuals, so varying an
#' endogenous regressor traces out the structural relationship.
#'
#' @param object An `"ivcf"` object.
#' @param newdata Optional data frame.
#' @param type `"link"` and `"response"` use the coefficients on the object's
#'   own normalisation; `"asf"` returns the average structural function
#'   whichever normalisation was used.
#' @param coefs Optional coefficient vector.
#' @param ... Ignored.
#' @export
predict.ivcf <- function(object, newdata = NULL,
                         type = c("link", "response", "asf"),
                         coefs = NULL, ...) {
  type <- match.arg(type)
  b <- if (is.null(coefs)) object$coefficients else coefs
  X <- aiv_design(object, newdata)

  is_v <- names(b) %in% colnames(object$residuals_first_stage)
  bs   <- b[!is_v]
  eta  <- as.vector(X[, names(bs), drop = FALSE] %*% bs)

  if (type == "asf" && object$normalize == "conditional") {
    lam <- b[is_v]
    eta <- eta / sqrt(1 + as.numeric(t(lam) %*% object$Sigma_v %*% lam))
  }
  if (type == "link") eta else object$family$linkinv(eta)
}

#' Wald test that the first-stage residuals have zero coefficients
#'
#' Under the model's assumptions this is a test of exogeneity of the
#' endogenous regressors. It requires no correction for estimation of the
#' first stage, because that correction is proportional to the coefficients
#' being tested.
#'
#' @param object An `"ivcf"` object.
#' @return A list with the Wald statistic, degrees of freedom and p-value.
#' @export
exogeneity_test <- function(object) {
  stopifnot(inherits(object, "ivcf"))
  i <- names(coef(object)) %in% colnames(object$residuals_first_stage)
  b <- coef(object)[i]
  V <- object$vcov[i, i, drop = FALSE]
  stat <- as.numeric(t(b) %*% solve(V, b))
  list(statistic = stat, df = length(b),
       p.value = stats::pchisq(stat, length(b), lower.tail = FALSE))
}

#' @export
print.ivcf <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nControl function (Rivers-Vuong), n = ", x$n,
      "\nNormalisation: ", x$normalize,
      if (x$normalize == "marginal")
        paste0("  (coefficients scaled by ", format(x$scale, digits = digits), ")"),
      "\n", sep = "")
  invisible(x)
}

#' @export
summary.ivcf <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z  <- object$coefficients / se
  structure(list(
    coefficients = cbind(Estimate = object$coefficients, `Std. Error` = se,
                         `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))),
    object = object), class = "summary.ivcf")
}

#' @export
print.summary.ivcf <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"),
                               ...) {
  o <- x$object
  cat("\nCall:\n", paste(deparse(o$call), collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, has.Pvalue = TRUE,
                      P.values = TRUE)

  item <- function(k, v) cat(formatC(k, width = 22, flag = "-"), v, "\n", sep = "")
  cat("\n")
  item("Estimator:", "control function (Rivers-Vuong)")
  item("Endogenous:", paste(o$endogenous, collapse = ", "))
  item("Instruments:", paste(o$instruments, collapse = ", "))
  item("Observations:", format(o$n, big.mark = ","))
  item("Normalisation:", if (o$normalize == "marginal")
    "Var(U) = 1 (marginal)" else "Var(U | Z, V) = 1 (conditional)")
  item("Scale factor:", format(o$scale, digits = digits))

  et <- exogeneity_test(o)
  item("Exogeneity test:", sprintf("chi2(%d) = %s, p = %s", et$df,
                                   format(et$statistic, digits = digits),
                                   format.pval(et$p.value, digits = digits)))
  invisible(x)
}
