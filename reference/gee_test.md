# Generalized Estimating Equations under technical replication with robust and non-robust Wald and Rao (score) tests

Generalized Estimating Equations under technical replication with robust
and non-robust Wald and Rao (score) tests

## Usage

``` r
gee_test(
  ...,
  use_geeasy = TRUE,
  use_jack_se = FALSE,
  cluster_corr_coef = NULL,
  skip_gee = FALSE
)
```

## Arguments

- ...:

  Arguments that you would pass to a regular
  [`geepack::geeglm`](https://rdrr.io/pkg/geepack/man/geeglm.html) call.
  Any observations with `NA` values in the data (response or covariates
  or id) will be dropped.

- use_geeasy:

  When TRUE, uses `geeasy` for gee estimation, when FALSE uses `geepack`

- use_jack_se:

  When TRUE uses jackknife standard errors (which take longer), when
  FALSE uses sandwich standard errors

- cluster_corr_coef:

  Optional within-cluster correlation coefficient. This will only be
  used when parameter estimation with a GEE fails and estimation must
  instead be performed with a GLM.

- skip_gee:

  When TRUE doesn't try to optimize with a GEE (just uses a GLM). This
  should only be used internally for testing.

## Examples

``` r
cars$id <- rep(1:5, each = 10)
gee_test(dist ~ speed, data = cars, family=poisson(link="log"), id = id)
#> Error in geelm(dist ~ speed, data = structure(list(speed = c(4, 4, 7,  : 
#>   could not find function "geelm"
#> Warning: GEE solver failed. Estimation will be done with a GLM, which will provide consistent parameter estimates but will not be as efficient as estimates from a GEE. Robust standard errors will be estimated using a clustered jackknife procedure.
#> Warning: Because estimation has been done with a GLM, there is no estimated cluster correlation coefficient. Therefore the robust score test cannot be run. In order to run a robust score test for this model, please input a within-cluster correlation coefficient to use.
#> 
#> Call:
#> gee_test(dist ~ speed, data = cars, family = poisson(link = "log"), 
#>     id = id)
#> 
#> 
#> Coefficient estimates:
#>                     Estimate Robust Std Error Lower 95% CI Upper 95% CI
#> (Intercept)       2.15096109       0.26098164    1.6394465    2.6624757
#> speed             0.09650242       0.01357674    0.0698925    0.1231123
#> correlation:alpha         NA               NA           NA           NA
#>                   Robust Wald p Robust Score p
#> (Intercept)        2.220446e-16             NA
#> speed              1.178058e-12             NA
#> correlation:alpha            NA             NA
```
