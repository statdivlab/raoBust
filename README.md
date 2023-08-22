
<!-- README.md is generated from README.Rmd. Please edit that file -->

# raoBust: Robust Rao (score) tests for Generalized Linear Models

`raoBust`, at its core, gives you all the important information from
`glm()`, but also with model misspecification-robust Rao tests (also
called score tests) and Wald tests.

Robust score tests have outstanding error rate performance in small
samples, **and** when your data is not drawn from a parametric family
(i.e., always). It’s shocking to me (Amy) how well they perform. They
have a reputation for being conservative in small samples, but I would
argue that this is a *very good thing*.

For now (because it’s what I need for my work) I only have
implementations for Poisson GLMs with log links. If you have another
specific case you’d like to request, please let me know at
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

The software was written by Amy D Willis, with input, insight and
inspiration from former lab member David S Clausen. All errors are Amy’s
fault.
