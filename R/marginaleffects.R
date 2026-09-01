# Extension methods for the marginaleffects package, which recognises a 0/1
# regressor as binary and reports the discrete 0 -> 1 contrast.

get_coef.aiv <- function(model, ...) model$coefficients

set_coef.aiv <- function(model, coefs, ...) {
  model$coefficients[names(coefs)] <- coefs
  model
}

get_vcov.aiv <- function(model, vcov = NULL, ...) {
  if (isFALSE(vcov)) return(NULL)
  if (is.matrix(vcov)) return(vcov)
  model$vcov
}

get_predict.aiv <- function(model, newdata = NULL, type = "response", ...) {
  est <- stats::predict(model, newdata = newdata,
                        type = if (identical(type, "link")) "link" else "response")
  data.frame(rowid = seq_along(est), estimate = as.numeric(est))
}

# Under the Rivers-Vuong normalisation Phi(x'beta) is P(y = 1 | X, V = 0), so
# the coefficients, their variance and the average structural function would
# sit on different scales. Refuse rather than rescale silently.
require_marginal <- function(model) {
  if (!identical(model$normalize, "marginal"))
    stop("marginal effects require normalize = \"marginal\".\n",
         "Under the Rivers-Vuong normalisation Phi(x'beta) is P(y = 1 | X, ",
         "V = 0), not the average structural function.\n",
         "Refit with normalize = \"marginal\", or use ",
         "predict(fit, type = \"asf\").", call. = FALSE)
}

get_coef.ivcf <- function(model, ...) {
  require_marginal(model)
  model$coefficients
}

set_coef.ivcf <- function(model, coefs, ...) {
  model$coefficients[names(coefs)] <- coefs
  model
}

get_vcov.ivcf <- function(model, vcov = NULL, ...) {
  if (isFALSE(vcov)) return(NULL)
  if (is.matrix(vcov)) return(vcov)
  require_marginal(model)
  model$vcov
}

get_predict.ivcf <- function(model, newdata = NULL, type = "response", ...) {
  require_marginal(model)
  est <- stats::predict(model, newdata = newdata,
                        type = if (identical(type, "link")) "link" else "response")
  data.frame(rowid = seq_along(est), estimate = as.numeric(est))
}

.onLoad <- function(libname, pkgname) {
  classes <- c("aiv", "ivcf")
  options(marginaleffects_model_classes =
            union(getOption("marginaleffects_model_classes", character()),
                  classes))

  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    ns <- asNamespace("marginaleffects")
    for (g in c("get_coef", "set_coef", "get_vcov", "get_predict"))
      for (cls in classes)
        registerS3method(g, cls, get(paste0(g, ".", cls)), envir = ns)
  }
  invisible()
}
