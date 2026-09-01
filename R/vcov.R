# Asymptotic variance of the profiled estimator.

# Asymptotic variance of the profiled estimator, evaluated at the fitted index
# with the auxiliary coefficients set to zero. Notation follows the paper:
# L_ab are blocks of the Hessian, G and H their partialled-out counterparts.
# B collects the coefficients of the two influence-function representations,
# so B S B' is the joint variance of (beta_hat, alpha_hat).
aiv_variance <- function(y, X_end, X_ex, Z_ex, family, beta, alpha, Omega) {
  n <- length(y)
  d <- index_derivs(family, y, as.vector(X_end %*% beta + X_ex %*% alpha))

  # Weighted cross-products: cp(A, B, w) = A' diag(w) B / n.
  cp <- function(A, B, w) crossprod(A * w, B) / n

  L_gg <- cp(Z_ex, Z_ex,  d$l2)
  L_gb <- cp(Z_ex, X_end, d$l2)
  L_ga <- cp(Z_ex, X_ex,  d$l2)
  L_aa <- cp(X_ex, X_ex,  d$l2)
  L_ab <- cp(X_ex, X_end, d$l2)

  S_gg <- cp(Z_ex, Z_ex, d$l1^2)
  S_aa <- cp(X_ex, X_ex, d$l1^2)
  S_ga <- cp(Z_ex, X_ex, d$l1^2)

  L_aa_inv <- solve(L_aa)
  A <- L_ga %*% L_aa_inv
  G <- L_gb - A %*% L_ab
  H <- L_gg - A %*% t(L_ga)
  W <- solve(H) %*% Omega %*% solve(H)

  M <- solve(t(G) %*% W %*% G, t(G) %*% W)
  P <- L_aa_inv %*% L_ab %*% M

  B <- rbind(cbind(-M, M %*% A),
             cbind( P, -(P %*% A + L_aa_inv)))
  S <- rbind(cbind(S_gg, S_ga),
             cbind(t(S_ga), S_aa))

  V <- B %*% S %*% t(B) / n
  dimnames(V) <- list(c(colnames(X_end), colnames(X_ex)),
                      c(colnames(X_end), colnames(X_ex)))

  # Smallest singular value of G is the generalised relevance condition: a
  # value near zero means a flat objective and imprecise estimates.
  list(vcov = V, relevance = min(svd(G)$d))
}

#' Test the relevance of the endogenous regressors
#'
#' Tests that the coefficients on the endogenous regressors are jointly zero —
#' the paper's test of regressor relevance. Note this concerns the effect of
#' the endogenous regressors on the OUTCOME, and is unrelated to instrument
#' relevance, which is the correlation between the instruments and the
#' endogenous regressors. The latter is reported separately by [aiv()] as
#' `relevance`, the smallest singular value of the partialled-out Jacobian.
#'
#' The variance is evaluated at the null rather than at the estimate. A Wald
#' statistic whose variance is evaluated at the estimate is subject to
#' the Hauck-Donner effect: it need not increase as the estimate moves further
#' from the null, so power can fall away for large effects. In a binary choice
#' model this happens once fitted probabilities approach zero and one, where
#' the estimated information collapses and the standard error is inflated.
#' Evaluating the variance at the null avoids it, as in `VGAM::wald.stat()`.
#'
#' Under the null the endogenous regressors drop out and, with the auxiliary
#' coefficients also zero, the model reduces to a generalised linear model of
#' the outcome on the exogenous regressors. `refit = TRUE` re-estimates the
#' remaining coefficients under that restriction, which is the recommended
#' default; `refit = FALSE` keeps them at their unrestricted values, which is
#' cheaper but evaluates the variance at a parameter point that is not the
#' restricted estimate.
#'
#' @param object An `"aiv"` object.
#' @param refit Re-estimate the exogenous coefficients under the null.
#' @return A list with the statistic, degrees of freedom, p-value, the
#'   standard errors used, and the unrestricted standard errors for comparison.
#' @references Yee, T. W. (2022). On the Hauck-Donner effect in Wald tests.
#'   \emph{Journal of the American Statistical Association} 117, 1763--1774.
#' @export
regressor_relevance_test <- function(object, refit = TRUE) {
  stopifnot(inherits(object, "aiv"))
  p   <- object$parse
  fam <- object$family
  b   <- object$beta_end

  alpha0 <- if (refit) {
    stats::glm.fit(p$X_ex, p$y, family = fam, intercept = FALSE)$coefficients
  } else {
    object$alpha
  }
  alpha0[is.na(alpha0)] <- 0

  V0 <- aiv_variance(p$y, p$X_end, p$X_ex, p$Z_ex, fam,
                     beta = rep(0, length(b)), alpha = alpha0,
                     Omega = crossprod(p$Z_ex) / object$n)
  Vb <- V0$vcov[names(b), names(b), drop = FALSE]

  stat <- drop(t(b) %*% solve(Vb, b))
  list(statistic = stat, df = length(b),
       p.value = stats::pchisq(stat, length(b), lower.tail = FALSE),
       se = sqrt(diag(Vb)),
       se_unrestricted = if (is.null(object$vcov)) NA_real_ else
         sqrt(diag(object$vcov)[names(b)]),
       refit = refit)
}
