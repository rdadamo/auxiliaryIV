#' @keywords internal
"_PACKAGE"

# Generics for which the package defines methods must be imported, otherwise
# the namespace cannot be loaded with only its stated dependencies.
#' @importFrom stats binomial coef family formula glm.control nobs predict
#' @importFrom stats terms vcov
NULL
