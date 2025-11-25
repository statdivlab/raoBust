# raoBust: Robust Rao (score) tests for Generalized Linear Models

`raoBust`, at its core, gives you all the important information from
[`glm()`](https://rdrr.io/r/stats/glm.html), but also with robust score
tests. Robust score tests are robust (in error rate control) to many
forms of model misspecification.

Robust score tests have outstanding error rate performance in small
samples, **and** when your data is not drawn from a parametric family
(i.e., always). It is shocking how well they perform. They are generally
conservative in small samples, which is a very good thing. You *should*
err on conservative when you have few samples. Most other tests are
anti-conservative in small samples.

We currently have score tests for coefficients in Poisson GLMs (log
link), Binomial GLMs (logit link), and Multinomial GLMs (log link),
including for linear combinations of parameters and simultaneous testing
(“ANOVA”). If you have another specific case you’d like to request,
please let us know at
[Issues](https://github.com/statdivlab/raoBust/issues) and label it as a
“feature request”.

This package also implements robust Wald tests. Many other packages
implement robust Wald tests – it is general methodology for robust score
tests that are the unique contribution of this package.

## Installation

You can install the development version of raoBust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("statdivlab/raoBust")
```

If it would be helpful for you for this package to be on CRAN, please
let us know.

## Example

This is a really silly example to show you the syntax. It looks just
like [`glm()`](https://rdrr.io/r/stats/glm.html) on the input side, but
outputs a summary table that includes robust Wald and Rao tests as well
as others.

``` r
library(raoBust)
#> Registered S3 methods overwritten by 'geeasy':
#>   method       from   
#>   drop1.geeglm MESS   
#>   drop1.geem   MESS   
#>   plot.geeglm  geepack
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

## Documentation

We additionally have a `pkgdown`
[website](https://statdivlab.github.io/raoBust/) that contains pre-built
versions of our function
[documentation](https://statdivlab.github.io/raoBust/reference/index.html)
and an introductory
[vignette](https://statdivlab.github.io/raoBust/articles/intro_raoBust.html).

## People

- Creator, maintainer: Amy D Willis
- Author: Sarah Teichman
- Author: David S Clausen
- Author: Shirley Mathur

All errors are Amy’s fault.

## Naming

Score tests were pioneered by C.R. Rao, an outstanding statistician and
methodologist. Score tests are sometimes, but increasingly rarely,
called Rao tests. While Rao’s work focused on non-robust score tests
(not robust score tests) it is in tribute to Rao that this package is
named. Many thanks to David Clausen for proposing and sharing the clever
portmanteau \`raoBust’ for a package implementing robust score tests.

**Comment from Amy** People of Color are consistently undervalued in
science and mathematics, including in statistical methodology. The
foundations of statistical methodology were built on the desire to
quantitatively show the inferiority of People of Color and Jews. I
believe that white supremacy continues to manifest in what names we give
methods: we are more likely to name a method after a person if they are
White, and more likely to name a method after its purpose, details or
another acronym if the developer was a Person of Color. For this reason,
I try to call score tests “Rao tests”. I invite you to join me.
