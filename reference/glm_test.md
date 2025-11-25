# Generalized Linear Models with robust and non-robust Wald and Rao (score) tests

Generalized Linear Models with robust and non-robust Wald and Rao
(score) tests

## Usage

``` r
glm_test(...)
```

## Arguments

- ...:

  Arguments that you would pass to a regular `glm` call. Any
  observations with `NA` values in the data (response or covariates)
  will be dropped.

## Examples

``` r
glm_test(dist ~ speed, data = cars, family=poisson(link="log"))
#> 
#> Call:
#> glm(dist ~ speed, data = cars, family = poisson(link = "log"))
#> 
#> 
#> Coefficient estimates:
#>               Estimate Non-robust Std Error Robust Std Error Lower 95% CI
#> (Intercept) 2.15096109          0.081774352      0.180014289     1.798140
#> speed       0.09650242          0.004404885      0.009234056     0.078404
#>             Upper 95% CI Non-robust Wald p Robust Wald p Robust Score p
#> (Intercept)    2.5037826     1.743527e-152             0   4.945013e-07
#> speed          0.1146008     2.177435e-106             0   3.956496e-05
```
