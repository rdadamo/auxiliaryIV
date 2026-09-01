# auxiv

Auxiliary instrumental variable (AIV) estimation for non-linear models with
endogenous regressors, implementing the estimator of D'Adamo, Weidner and
Windmeijer (2026).

The excluded instruments are added to the model as auxiliary regressors, and
the coefficients on the endogenous regressors are chosen so that the estimated
auxiliary coefficients come out as close to zero as possible. Unlike the
control function or bivariate probit, the estimator imposes no functional form
or distributional assumption on the first stage. It is consistent when the
coefficients on the endogenous regressors are zero, and locally sign consistent
elsewhere.

## Installation

```r
# install.packages("remotes")
remotes::install_github("rdadamo/auxiliaryIV")
```

## Usage

The formula follows the `ivreg` convention: regressors on the left of the `|`,
the complete instrument set on the right — including any exogenous regressors
that appear in both.

```r
library(auxiv)

set.seed(1)
n <- 500
z <- rnorm(n)                                  # instrument
w <- rnorm(n)                                  # exogenous regressor
u <- rnorm(n)                                  # structural error
x <- 0.8 * z + 0.4 * w + 0.7 * u + rnorm(n)    # endogenous regressor
y <- as.numeric(0.3 - 0.4 * w + u >= 0)        # outcome, x has no effect

fit <- aiv(y ~ x + w | z + w, data = data.frame(y, x, w, z))
summary(fit)
```

```
Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  0.23137    0.05883   3.933 8.38e-05 ***
x            0.04089    0.06773   0.604    0.546
w           -0.49288    0.06560  -7.513 5.77e-14 ***

Family:               binomial (probit)
Endogenous:           x
Instruments:          z
Observations:         500
Identification:       exactly identified (1 instrument, 1 endogenous regressor)
Objective:            3.767e-18
Auxiliary coef.:      z = 1.92e-09
Instrument relevance: 0.4921

Endogenous regressor relevance test
  Null:               x = 0
  Statistic:          chi2(1) = 0.3648, p = 0.5459
  Variance:           evaluated at the null

Numerical optimisation
  Optimiser:          grid search + optimize
  Search range:       [-2.533, 2.533]
  Model fits on:      [-1.055, 2.322], 17 of 19 points
  Local minima:       1
```

## What the package provides

| Function | Purpose |
|---|---|
| `aiv()` | The AIV estimator. Probit, logit, cloglog and Poisson through the standard `family` interface. Standard errors valid under over-identification. |
| `regressor_relevance_test()` | Tests that the endogenous regressors have no effect on the outcome, with the variance evaluated at the null. Reported by `summary()`. |
| `ivcf()` | Rivers–Vuong control function, for comparison, under either normalisation of the structural error variance. |
| `exogeneity_test()` | Tests that the first-stage residuals have zero coefficients in the control function fit. |
| `aiv_control()` | Settings for the search over the endogenous coefficients. |

Average partial effects are delegated to
[marginaleffects](https://marginaleffects.com), which recognises a 0/1
regressor and reports the discrete contrast rather than a derivative:

```r
marginaleffects::avg_slopes(fit, variables = "x")
```

## Citation

```r
citation("auxiv")
```

## Status

Under development, alongside the paper. The interface may still change.

## Licence

GPL-3
