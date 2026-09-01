# Auxiliary IV estimator: parsing, estimation, and the main entry point.

# Parse a two-part formula in the ivreg convention,
# y ~ x_end + x_ex | z_ex + x_ex, whose second part lists the complete
# instrument set. Endogenous regressors are those among the regressors but not
# among the instruments.
aiv_parse <- function(formula, data) {
  f <- Formula::as.Formula(formula)
  if (length(f)[2] != 2L)
    stop("'formula' needs two right-hand parts: y ~ regressors | instruments")

  mf <- stats::model.frame(f, data = data)
  y  <- stats::model.response(mf)
  if (is.factor(y)) y <- as.numeric(y) - 1

  X <- stats::model.matrix(f, data = mf, rhs = 1)
  Z <- stats::model.matrix(f, data = mf, rhs = 2)

  end <- setdiff(colnames(X), colnames(Z))
  ex  <- intersect(colnames(X), colnames(Z))
  iv  <- setdiff(colnames(Z), colnames(X))

  if (length(end) == 0L)
    stop("no endogenous regressors: every regressor is also an instrument")
  if (length(iv) < length(end))
    stop("order condition fails: ", length(iv), " excluded instrument(s) for ",
         length(end), " endogenous regressor(s)")

  list(y = as.numeric(y), X = X,
       X_end = X[, end, drop = FALSE],
       X_ex  = X[, ex,  drop = FALSE],
       Z_ex  = Z[, iv,  drop = FALSE],
       end = end, ex = ex, iv = iv,
       terms = stats::terms(mf),
       xlevels = stats::.getXlevels(stats::terms(mf), mf),
       model_frame = mf)
}

# Design matrix for prediction: the fitted one, or rebuilt from newdata.
aiv_design <- function(object, newdata) {
  if (is.null(newdata)) return(object$parse$X)
  tt <- stats::delete.response(object$terms)
  stats::model.matrix(tt, stats::model.frame(tt, newdata,
                                             xlev = object$parse$xlevels))
}

# Outer objective: for given endogenous coefficients, fit the model with the
# instruments added as auxiliary regressors and return the weighted squared
# norm of their estimated coefficients. Values of beta large enough to separate 
# the data are not admissible and return Inf.
aiv_objective <- function(beta, y, X_end, D, k_iv, family, Omega, glm_control) {
  fit <- suppressWarnings(stats::glm.fit(
    x = D, y = y, family = family, offset = as.vector(X_end %*% beta),
    intercept = FALSE, control = glm_control))

  gamma <- fit$coefficients[seq_len(k_iv)]
  if (!fit$converged || anyNA(gamma)) return(Inf)
  drop(crossprod(gamma, Omega %*% gamma))
}

# Minimises the objective over beta. One endogenous regressor: grid search
# from zero, refined by optimize(). Multiple endogenous regressors: the 
# optimiser named by optfct, from start.
aiv_outer <- function(obj, X_end, control, start = NULL) {
  sdx <- apply(X_end, 2, stats::sd)
  sdx[!is.finite(sdx) | sdx <= 0] <- 1
  half <- control$range / sdx

  if (ncol(X_end) > 1L) {
    if (is.null(start)) start <- rep(0, ncol(X_end))
    if (control$optfct == "nlminb") {
      r    <- stats::nlminb(start, obj)
      par  <- r$par
      value <- r$objective
      info <- list(optimiser = "nlminb", convergence = r$convergence,
                   counts = r$evaluations, message = r$message)
    } else {
      r    <- stats::optim(start, obj, method = "Nelder-Mead",
                           control = list(reltol = control$tol))
      par  <- r$par
      value <- r$value
      info <- list(optimiser = "optim (Nelder-Mead)",
                   convergence = r$convergence, counts = r$counts,
                   message = r$message)
    }
    if (!is.finite(value))
      stop("the model could not be fitted at any candidate value")
    return(list(beta = par, objective = value,
                objective_zero = obj(rep(0, ncol(X_end))),
                algoInfo = info, grid = NULL, values = NULL, minima = NULL))
  }

  # Beta is only fittable in an interval around zero — further out the index
  # separates the data. Scan outward, stopping each side at the first failure.
  m    <- (control$n_grid - 1) %/% 2
  step <- half / m
  grid <- 0
  val  <- obj(0)
  for (side in c(1, -1)) {
    for (k in seq_len(m)) {
      b <- side * k * step
      v <- obj(b)
      grid <- c(grid, b)
      val  <- c(val, v)
      if (!is.finite(v)) break
    }
  }
  o <- order(grid); grid <- grid[o]; val <- val[o]

  if (all(!is.finite(val)))
    stop("the model could not be fitted at any candidate value")

  # Find every local minimum on the grid, i.e. each point whose objective is
  # lower than at the two points either side of it.
  vv <- val
  vv[!is.finite(vv)] <- Inf
  cand <- which(is.finite(val) &
                vv <  c(Inf, vv[-length(vv)]) &
                vv <= c(vv[-1L], Inf))
  if (length(cand) == 0L) cand <- which.min(vv)

  # Improves one grid minimum by searching between the two grid points either
  # side of it. The interval has to span a point where the model did not fit,
  # since the minimum can lie between the last one that fitted and the first
  # that did not. optimize() handles the infinite objective there but warns
  # about it, which is not something the user can act on, so it is silenced.
  refine <- function(i) {
    lo <- grid[max(1L, i - 1L)]
    hi <- grid[min(length(grid), i + 1L)]
    if (hi <= lo) return(c(grid[i], val[i]))
    o <- suppressWarnings(
      stats::optimize(obj, interval = c(lo, hi), tol = control$tol))
    if (is.finite(o$objective) && o$objective < val[i]) c(o$minimum, o$objective)
    else c(grid[i], val[i])
  }
  sol    <- vapply(cand, refine, numeric(2))
  minima <- data.frame(beta = sol[1, ], objective = sol[2, ])

  # Report the global minimum. Where several are indistinguishable the
  # objective cannot choose between them, so fall back on the one closest to
  # zero.
  # Indistinguishable is judged against the objective at beta = 0, the same
  # yardstick the no-root check in aiv() uses.
  obj0 <- val[grid == 0]
  tied <- which(minima$objective - min(minima$objective) <=
                  1e-6 * max(obj0, .Machine$double.eps))
  k    <- if (length(tied) > 1L) tied[which.min(abs(minima$beta[tied]))] else tied

  if (length(tied) > 1L)
    warning(length(tied), " local minima give essentially the same objective, ",
            "at beta = ",
            paste(format(minima$beta[tied], digits = 4), collapse = ", "),
            ", so the estimate is not pinned down. Reporting the one closest ",
            "to zero; see the 'minima' element for all of them.",
            call. = FALSE)

  # The search can end in three ways: at an interior minimum, at the edge of
  # the search range, or at the edge of the region where the model can be
  # fitted. The second and third cases trigger warnings.
  i    <- cand[k]
  next_to_unfitted <- !is.finite(val[max(1L, i - 1L)]) ||
                      !is.finite(val[min(length(val), i + 1L)])

  if (abs(grid[i]) >= half * (1 - 1e-8))
    warning("the estimate is at the edge of the search range (beta = ",
            format(minima$beta[k], digits = 4), ", range = ",
            format(control$range, digits = 4),
            "). Widen it with control = aiv_control(range = ...).",
            call. = FALSE)
  else if (next_to_unfitted)
    warning("the estimate is at the edge of the region where the model can be ",
            "fitted (beta = ", format(minima$beta[k], digits = 4),
            "). This is a boundary solution rather than an interior minimum, ",
            "so the usual standard errors may not apply.", call. = FALSE)

  list(beta = minima$beta[k], objective = minima$objective[k],
       objective_zero = val[grid == 0],
       algoInfo = list(optimiser = "grid search + optimize",
                       range = control$range, half = half,
                       counts = length(val),
                       admissible = sum(is.finite(val)),
                       admissible_range = range(grid[is.finite(val)]),
                       local_minima = nrow(minima)),
       grid = grid, values = val, minima = minima)
}

#' Control parameters for aiv()
#'
#' @param range Half-width of the grid, in standard deviations of the
#'   endogenous regressor. A coefficient outside this range would move the
#'   index by more than `range` standard deviations, which no binary choice
#'   model can support. Used only with a single endogenous regressor.
#' @param n_grid Number of grid points. The grid only has to bracket the
#'   minimum; `optimize()` refines within one step of it. Single endogenous
#'   regressor only.
#' @param tol Tolerance of the outer optimiser.
#' @param optfct Outer optimiser used when there is more than one endogenous
#'   regressor, following the convention of the `gmm` package. `"nlminb"` is
#'   a quasi-Newton method and is the default; `"optim"` is Nelder-Mead, which
#'   is more robust to a badly behaved objective but stalls well short of the
#'   optimum once its simplex collapses. With a single endogenous regressor
#'   the grid search is used instead and this is ignored.
#' @param glm A list from [stats::glm.control()] governing the inner fit.
#' @return A list of control parameters.
#' @export
aiv_control <- function(range = 4, n_grid = 25, tol = 1e-6,
                        optfct = c("nlminb", "optim"),
                        glm = glm.control()) {
  list(range = range, n_grid = n_grid, tol = tol,
       optfct = match.arg(optfct), glm = glm)
}

#' Auxiliary instrumental variable estimation
#'
#' Estimates a non-linear model with endogenous regressors using the auxiliary
#' instrumental variable (AIV) estimator of D'Adamo, Weidner and Windmeijer
#' (2026).
#'
#' The excluded instruments are added to the model as auxiliary regressors and
#' the coefficients on the endogenous regressors are chosen so that the
#' estimated auxiliary coefficients are as close to zero as possible.
#'
#' The auxiliary coefficients and the coefficients on the exogenous regressors
#' are concentrated out in a concave inner problem, leaving a search over the
#' endogenous coefficients alone.
#'
#' The estimator imposes no functional form or distributional assumption on the
#' first stage. It is consistent when the coefficients on the endogenous
#' regressors are zero, and locally sign consistent elsewhere.
#'
#' @param formula A two-part formula, `y ~ x_end + x_ex | z_ex + x_ex`, whose
#'   second part lists the complete instrument set.
#' @param data A data frame.
#' @param family A `family` object; the binomial links and `poisson()` are
#'   supported.
#' @param start Optional starting value, used only when there is more than one
#'   endogenous regressor.
#' @param control A list from [aiv_control()], or a plain list of the settings
#'   to change, such as `list(range = 8)`.
#' @param se Whether to compute the asymptotic variance, and hence the
#'   standard errors.
#' @param Omega Weight matrix for the auxiliary coefficients. Defaults to
#'   \eqn{Z'Z/n} for the excluded instruments.
#' @return An object of class `"aiv"`. `gamma` holds the auxiliary
#'   coefficients at the estimate, the quantity being driven to zero. The
#'   profile of the outer objective is returned in `grid` and `grid_values`,
#'   and every local minimum found in `minima`; where there is more than one
#'   the estimate is not uniquely defined and the one closest to zero is
#'   reported, with a warning.
#' @references
#' D'Adamo, R., M. Weidner and F. Windmeijer (2026). Auxiliary IV estimation
#' for nonlinear models. Working paper.
#' @seealso [regressor_relevance_test()] for the test that the endogenous
#'   regressors have no effect, [ivcf()] for the Rivers-Vuong control function
#'   estimator, and [aiv_control()] for the settings of the search.
#' @examples
#' set.seed(1)
#' n <- 500
#' z <- rnorm(n); # instrument
#' w <- rnorm(n); # exogenous regressor
#' u <- rnorm(n); # structural error
#' x <- 0.8 * z + 0.4 * w + 0.7 * u + rnorm(n) # endogenous regressor
#' y <- as.numeric(0.3 - 0.4 * w + u >= 0) # outcome
#' # Estimate model with auxiliary IV estimator
#' fit <- aiv(y ~ x + w | z + w, data = data.frame(y, x, w, z))
#' summary(fit)
#' @export
aiv <- function(formula, data, family = binomial("probit"),
                start = NULL, control = aiv_control(), se = TRUE,
                Omega = NULL) {
  cl <- match.call()
  if (is.character(family)) family <- get(family, mode = "function")()
  if (is.function(family))  family <- family()

  # As glm() does, fill a plain list out with the defaults, so that
  # control = list(range = 8) works without calling aiv_control().
  control <- do.call(aiv_control, control)

  p <- aiv_parse(formula, data)
  n <- length(p$y)
  k_iv <- ncol(p$Z_ex)

  D <- cbind(p$Z_ex, p$X_ex)
  if (is.null(Omega)) Omega <- crossprod(p$Z_ex) / n
  stopifnot(nrow(Omega) == k_iv, ncol(Omega) == k_iv)

  obj <- function(b)
    aiv_objective(b, p$y, p$X_end, D, k_iv, family, Omega, control$glm)

  outer_fit <- aiv_outer(obj, p$X_end, control, start)
  beta <- outer_fit$beta
  names(beta) <- p$end

  # Auxiliary coefficients at the estimate: the quantity being driven to zero.
  # Exactly zero only under just identification; over-identified, a positive
  # objective is the norm and its size is a specification test.
  fit_aux <- stats::glm.fit(x = D, y = p$y, family = family,
                            offset = as.vector(p$X_end %*% beta),
                            intercept = FALSE, control = control$glm)
  gamma <- fit_aux$coefficients[seq_len(k_iv)]
  names(gamma) <- p$iv

  # Just identified, gamma = 0 has as many equations as unknowns and a root
  # should exist. An objective that stays well above its value at zero means
  # the moment condition cannot be satisfied anywhere the model can be fitted,
  # so the estimate minimises a criterion rather than solving it. The relative
  # comparison is against beta = 0 because the objective carries the scale of
  # Omega; a genuine numerical root falls many orders of magnitude below it.
  # This needs no grid, so it applies with any number of endogenous regressors.
  obj0 <- outer_fit$objective_zero
  if (k_iv == ncol(p$X_end) && length(obj0) == 1L && is.finite(obj0) &&
      obj0 > 0 && outer_fit$objective > 1e-6 * obj0)
    warning("the auxiliary coefficients cannot be driven to zero anywhere ",
            "in the admissible range: objective ",
            format(outer_fit$objective, digits = 3), " against ",
            format(obj0, digits = 3), " at beta = 0, with as many ",
            "instruments as endogenous regressors. The estimate minimises ",
            "the objective rather than solving the moment condition.",
            call. = FALSE)

  # Coefficients on the exogenous regressors, with the auxiliary coefficients
  # set to zero.
  fit_ex <- stats::glm.fit(x = p$X_ex, y = p$y, family = family,
                           offset = as.vector(p$X_end %*% beta),
                           intercept = FALSE, control = control$glm)
  alpha <- fit_ex$coefficients
  alpha[is.na(alpha)] <- 0
  names(alpha) <- p$ex

  V <- if (se)
    aiv_variance(p$y, p$X_end, p$X_ex, p$Z_ex, family, beta, alpha, Omega)

  coefs <- c(beta, alpha)[colnames(p$X)]

  structure(list(
    coefficients = coefs,
    beta_end = beta, alpha = alpha,
    vcov = if (se) V$vcov[names(coefs), names(coefs), drop = FALSE],
    relevance = if (se) V$relevance,
    gamma = gamma,
    objective = outer_fit$objective,
    objective_zero = obj0,
    algoInfo = outer_fit$algoInfo,
    grid = outer_fit$grid, grid_values = outer_fit$values,
    minima = outer_fit$minima,
    family = family, formula = formula, terms = p$terms, call = cl,
    endogenous = p$end, exogenous = p$ex, instruments = p$iv,
    n = n, k_end = ncol(p$X_end), k_iv = k_iv,
    overidentified = k_iv > ncol(p$X_end),
    model = p$model_frame, parse = p
  ), class = "aiv")
}
