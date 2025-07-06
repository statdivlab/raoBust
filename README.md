
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raoBust: Robust Rao (score) tests for Generalized Linear Models

<!-- badges: start -->

[![R-CMD-check](https://github.com/statdivlab/raoBust/workflows/R-CMD-check/badge.svg)](https://github.com/statdivlab/raoBust/actions)
[![codecov](https://codecov.io/github/statdivlab/raoBust/coverage.svg?branch=main)](https://app.codecov.io/github/statdivlab/raoBust)
<!-- badges: end -->

`raoBust`, at its core, gives you all the important information from
`glm()`, but also with model misspecification-robust Rao tests (also
called score tests) and Wald tests.

Robust score tests have outstanding error rate performance in small
samples, **and** when your data is not drawn from a parametric family
(i.e., always). It is shocking how well they perform. They are
generally conservative in small samples, which is a very good thing. 
You *should* err on conservative when you have few samples. Most 
other tests are anti-conservative in small samples. 

We currently have Rao tests for coefficients in Poisson GLMs (log 
link), Binomial GLMs (logit link), and Multinomial GLMs (log link), 
including for linear combinations of parameters and simultaneous testing
("ANOVA"). If you have another
specific case you’d like to request, please let us know at
[Issues](https://github.com/statdivlab/raoBust/issues) and label it as a
“feature request”.

## Installation

You can install the development version of raoBust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("statdivlab/raoBust")
```

## Example

This is a really silly example to show you the syntax. It looks just
like `glm()` on the input side, but outputs a summary table that
includes robust Wald and Rao tests as well as others.

``` r
library(raoBust)
glm_test(dist ~ speed, data = cars, family=poisson(link="log"))
#>               Estimate Non-robust Std Error Robust Std Error Non-robust Wald p
#> (Intercept) 2.15096109          0.081774352      0.180014289     1.743527e-152
#> speed       0.09650242          0.004404885      0.009234056     2.177435e-106
#>             Robust Wald p Robust Score p
#> (Intercept)             0   0.0406316810
#> speed                   0   0.0000472766
```

## People

- Creator, maintainer: Amy D Willis
- Author: Sarah Teichman
- Author: David S Clausen
- Author: Shirley Mathur

All errors are Amy’s fault.
