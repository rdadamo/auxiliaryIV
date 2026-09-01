# Methods for fitted "aiv" objects.

#' @export
coef.aiv <- function(object, ...) object$coefficients

#' @export
vcov.aiv <- function(object, ...) object$vcov

#' @export
nobs.aiv <- function(object, ...) object$n

#' @export
formula.aiv <- function(x, ...) x$formula

#' @export
terms.aiv <- function(x, ...) x$terms

#' @export
family.aiv <- function(object, ...) object$family

#' Predictions from an AIV fit
#'
#' `coefs` allows prediction from a coefficient vector other than the fitted
#' one, which is what lets other packages obtain delta-method standard errors
#' by perturbing the coefficients and re-predicting.
#'
#' @param object An `"aiv"` object.
#' @param newdata Optional data frame.
#' @param type `"link"` for the index, `"response"` for the mean.
#' @param coefs Optional coefficient vector.
#' @param ... Ignored.
#' @export
predict.aiv <- function(object, newdata = NULL,
                        type = c("link", "response"), coefs = NULL, ...) {
  type <- match.arg(type)
  b <- if (is.null(coefs)) object$coefficients else coefs
  X <- aiv_design(object, newdata)
  eta <- as.vector(X[, names(b), drop = FALSE] %*% b)
  if (type == "link") eta else object$family$linkinv(eta)
}

aiv_ident <- function(x) {
  paste0(if (x$overidentified) "over-identified" else "exactly identified",
         " (", x$k_iv, " instrument", if (x$k_iv != 1L) "s", ", ",
         x$k_end, " endogenous regressor", if (x$k_end != 1L) "s", ")")
}

#' @export
print.aiv <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n", x$family$family, " (", x$family$link, "),  n = ", x$n, "\n",
      aiv_ident(x), "\n", sep = "")
  invisible(x)
}

#' @export
summary.aiv <- function(object, ...) {
  se <- if (is.null(object$vcov)) rep(NA_real_, length(object$coefficients))
        else sqrt(diag(object$vcov))
  z <- object$coefficients / se
  structure(list(
    coefficients = cbind(Estimate = object$coefficients, `Std. Error` = se,
                         `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z))),
    object = object), class = "summary.aiv")
}

#' @export
print.summary.aiv <- function(x, digits = max(3L, getOption("digits") - 3L),
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
  item("Family:", paste0(o$family$family, " (", o$family$link, ")"))
  item("Endogenous:", paste(o$endogenous, collapse = ", "))
  item("Instruments:", paste(o$instruments, collapse = ", "))
  item("Observations:", format(o$n, big.mark = ","))
  item("Identification:", aiv_ident(o))
  item("Objective:", format(o$objective, digits = digits))
  item("Auxiliary coef.:", paste(sprintf("%s = %s", names(o$gamma),
                                         format(o$gamma, digits = digits)),
                                 collapse = ", "))
  if (!is.null(o$relevance))
    item("Instrument relevance:", format(o$relevance, digits = digits))

  if (!is.null(o$vcov)) {
    rt <- try(regressor_relevance_test(o), silent = TRUE)
    if (!inherits(rt, "try-error")) {
      cat("\nEndogenous regressor relevance test\n")
      item("  Null:", paste0(paste(o$endogenous, collapse = " = "), " = 0"))
      item("  Statistic:", sprintf("chi2(%d) = %s, p = %s", rt$df,
                                   format(rt$statistic, digits = digits),
                                   format.pval(rt$p.value, digits = digits)))
      item("  Variance:", "evaluated at the null")
      # A large gap between the two standard errors indicates a Hauck-Donner
      # effect, i.e. an unrestricted Wald test with understated significance.
      r <- max(rt$se_unrestricted / rt$se)
      if (is.finite(r) && r > 1.5)
        cat("\nNote: the unrestricted standard error is ",
            format(r, digits = 2), " times the one evaluated at the null, ",
            "which\nindicates a Hauck-Donner effect; prefer the test above.\n",
            sep = "")
    }
  }

  # Reported rather than acted on, as in the gmm package: the user is told how
  # the outer problem was solved and can re-run with a different optimiser.
  a <- o$algoInfo
  if (!is.null(a)) {
    num <- function(v) format(v, digits = digits)
    cat("\nNumerical optimisation\n")
    item("  Optimiser:", a$optimiser)
    if (!is.null(a$half))
      item("  Search range:", sprintf("[%s, %s]", num(-a$half), num(a$half)))
    if (!is.null(a$admissible_range))
      item("  Model fits on:", sprintf("[%s, %s], %d of %d points",
                                       num(a$admissible_range[1]),
                                       num(a$admissible_range[2]),
                                       a$admissible, a$counts))
    else if (!is.null(a$counts)) {
      item("  Function eval.:", a$counts[1])
      if (length(a$counts) > 1L && !is.na(a$counts[2]))
        item("  Gradient eval.:", a$counts[2])
    }
    if (!is.null(a$convergence)) item("  Convergence code:", a$convergence)
    if (!is.null(a$message)) item("  Message:", a$message)
    if (!is.null(a$local_minima))
      item("  Local minima:", if (a$local_minima > 1L)
        paste0(a$local_minima, " at beta = ",
               paste(format(o$minima$beta, digits = digits), collapse = ", "))
        else a$local_minima)
  }
  invisible(x)
}
