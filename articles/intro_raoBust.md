# Introduction to raoBust

We’ll start by installing `raoBust` if we haven’t already, and then
loading it.

``` r
# if (!require("remotes", quietly = TRUE))
#     install.packages("remotes")
# 
# remotes::install_github("statdivlab/raoBust")

library(raoBust)
#> Registered S3 methods overwritten by 'geeasy':
#>   method       from   
#>   drop1.geeglm MESS   
#>   drop1.geem   MESS   
#>   plot.geeglm  geepack
library(geeasy) # also load geeasy, which raoBust uses in the back-end
#> Loading required package: geepack
```

## Introduction

`raoBust` is an R package that implements robust Wald tests and Rao
(score) tests for generalized linear models. We like robust tests
because they are robust (retain error rate control) to many forms of
model misspecification. We especially like robust Rao tests because they
have strong error rate control in small samples.

`raoBust` implements robust tests for coefficients in Poisson GLMs (log
link), Binomial GLMs (logit link), linear models (identity link), and
Multinomial GLMs (log link).

## Model fitting

### GLMs

We will demonstrate `raoBust` using the `mtcars` dataset.

``` r
head(mtcars)
#>                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
#> Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#> Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
```

First, let’s consider fitting a Poisson GLM to estimate and test the
expected fold-difference in `mpg` associated with a one-unit increase in
horsepower, account for type of transmission. We choose Poisson
regression because we want to estimate a fold-difference in means, and
because we are using robust tests we will not worry about whether we
actually expect `mpg` measurements to follow a Poisson distribution.

``` r
glm_test(formula = mpg ~ hp + am, data = mtcars, family = poisson(link = "log"))
#> 
#> Call:
#> glm(formula = mpg ~ hp + am, data = mtcars, family = poisson(link = "log"))
#> 
#> 
#> Coefficient estimates:
#>                 Estimate Non-robust Std Error Robust Std Error Lower 95% CI
#> (Intercept)  3.327082221         0.1151938176     0.0557891828  3.217737432
#> hp          -0.003110758         0.0006656492     0.0003679554 -0.003831938
#> am           0.234033015         0.0835981015     0.0479722881  0.140009058
#>             Upper 95% CI Non-robust Wald p Robust Wald p Robust Score p
#> (Intercept)  3.436427010     1.982292e-183  0.000000e+00   6.558146e-06
#> hp          -0.002389579      2.964427e-06  0.000000e+00   2.397949e-03
#> am           0.328056971      5.118157e-03  1.068933e-06   9.288962e-04
```

Looking at the output, we can see our estimated fold-difference in `mpg`
associated with a one-unit increase in horsepower is
`exp(-0.003) = 0.997`, with a robust Rao test p-value of $0.002$. We can
compare this to the results we would get if we ran a typical `glm`.

``` r
summary(glm(formula = mpg ~ hp + am, data = mtcars, family = poisson(link = "log")))
#> 
#> Call:
#> glm(formula = mpg ~ hp + am, family = poisson(link = "log"), 
#>     data = mtcars)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.3270822  0.1151938  28.882  < 2e-16 ***
#> hp          -0.0031108  0.0006656  -4.673 2.96e-06 ***
#> am           0.2340330  0.0835981   2.800  0.00512 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 54.524  on 31  degrees of freedom
#> Residual deviance: 10.727  on 29  degrees of freedom
#> AIC: Inf
#> 
#> Number of Fisher Scoring iterations: 4
```

Here we see the same estimates, but `glm` does not give us robust
standard errors or robust Wald and Rao tests.

### GEEs, accounting for clusters

In `raoBust` we can also use a GEE framework to account for clustered
observations. Let’s add a cluster variable to `mtcars`.

``` r
mtcars$cluster <- rep(1:8, 4)
gee_test(formula = mpg ~ hp + am, data = mtcars, family = poisson(link = "log"), id = cluster)
#> 
#> Call:
#> gee_test(formula = mpg ~ hp + am, data = mtcars, family = poisson(link = "log"), 
#>     id = cluster)
#> 
#> 
#> Coefficient estimates:
#>                      Estimate Robust Std Error Lower 95% CI Upper 95% CI
#> (Intercept)        3.32490495     0.0660228787  3.195502488  3.454307417
#> hp                -0.00309798     0.0003224075 -0.003729887 -0.002466073
#> am                 0.23472057     0.0560260439  0.124911541  0.344529597
#> correlation:alpha  0.09206717               NA           NA           NA
#>                   Robust Wald p Robust Score p
#> (Intercept)        0.000000e+00    0.005722124
#> hp                 0.000000e+00    0.016360642
#> am                 2.795817e-05    0.022199182
#> correlation:alpha            NA             NA
```

Again we have the same estimates, but different standard errors and
robust Wald and Rao test p-values once we account for clustering.

### Multinomial models

In `raoBust`, we have a function
[`multinom_test()`](../reference/multinom_test.md) to run robust tests
for multinomial regression. Important arguments are `strong` and `j`.
When we set `strong = TRUE`, we test the strong null hypothesis that
$\beta_{kj} = 0$ for all $k$ and $j$, where $k$ represents a covariate
and $j$ represents a level of the outcome variable. This means that
we’re testing the hypothesis that all non-intercept coefficients are
equal to $0$. When `strong = FALSE` then we must specify a outcome level
to test with the `j` argument. In this case, we test the hypothesis that
$\beta_{kj} = 0$ for the `j` specified and all covariates $k$.

`multinom_test` takes in a `Y` matrix, with `n` rows and `J` columns for
a dataset with `n` samples and a multi-category outcome with `J`
categories, and a design matrix `X` that includes metadata. We will
simulate data to demonstrate
[`multinom_test()`](../reference/multinom_test.md) on.

``` r
set.seed(123)  # for reproducibility

# design matrix: 30 samples x 2 covariates
n <- 30
p <- 2
J <- 5
X <- cbind(1, rep(0:1, each = n / 2), rnorm(n))

# coefficient matrix B: 3 rows x 5 columns
B <- matrix(rnorm(3 * J, mean = 0, sd = 0.5), nrow = 3, ncol = 5)
B[1, ] <- B[1, ] + 3

# generate Poisson outcomes
Y <- matrix(NA, nrow = n, ncol = J)
for (j in 1:J) {
  lambda <- exp(X %*% B[, j])
  Y[, j] <- rpois(n, lambda)
}
```

We can now fit our model. We’ll start by testing the strong null
hypothesis that all non-intercept coefficients are equal to $0$.

``` r
strong_test <- multinom_test(Y = Y, X = X, strong = TRUE)
strong_test$coef_tab
#>    Category   Covariate    Estimate Robust Std Error Lower 95% CI Upper 95% CI
#> 11        4           1 -1.52159024       0.14168311   -1.7992840  -1.24389645
#> 2         1           1 -1.19778008       0.10005117   -1.3938768  -1.00168340
#> 8         3           1 -1.13105847       0.11182994   -1.3502411  -0.91187582
#> 4         2 (intercept)  1.04369722       0.08435207    0.8783702   1.20902425
#> 7         3 (intercept)  0.93283338       0.08293878    0.7702764   1.09539040
#> 1         1 (intercept)  0.78509417       0.08361962    0.6212027   0.94898561
#> 5         2           1 -0.64777461       0.09662728   -0.8371606  -0.45838862
#> 9         3           2 -0.63990526       0.05227395   -0.7423603  -0.53745021
#> 12        4           2 -0.60215382       0.06820800   -0.7358390  -0.46846859
#> 10        4 (intercept)  0.47105039       0.12368367    0.2286349   0.71346592
#> 6         2           2 -0.13838624       0.03988563   -0.2165606  -0.06021185
#> 3         1           2 -0.04985085       0.04445023   -0.1369717   0.03727000
#>    Robust Wald p Robust Score p
#> 11  6.648178e-27             NA
#> 2   5.000849e-33             NA
#> 8   4.784250e-24             NA
#> 4   3.654085e-35             NA
#> 7   2.389181e-29             NA
#> 1   6.064411e-21             NA
#> 5   2.030015e-11             NA
#> 9   1.868292e-34             NA
#> 12  1.063748e-18             NA
#> 10  1.398072e-04             NA
#> 6   5.212790e-04             NA
#> 3   2.620759e-01             NA
strong_test$test_stat
#> [1] 25.10367
strong_test$p
#> [1] 0.001492897

weak_test <- multinom_test(Y = Y, X = X, strong = FALSE, j = 3)
weak_test$coef_tab
#>    Category   Covariate    Estimate Robust Std Error Lower 95% CI Upper 95% CI
#> 11        4           1 -1.52159024       0.14168311   -1.7992840  -1.24389645
#> 2         1           1 -1.19778008       0.10005117   -1.3938768  -1.00168340
#> 8         3           1 -1.13105847       0.11182994   -1.3502411  -0.91187582
#> 4         2 (intercept)  1.04369722       0.08435207    0.8783702   1.20902425
#> 7         3 (intercept)  0.93283338       0.08293878    0.7702764   1.09539040
#> 1         1 (intercept)  0.78509417       0.08361962    0.6212027   0.94898561
#> 5         2           1 -0.64777461       0.09662728   -0.8371606  -0.45838862
#> 9         3           2 -0.63990526       0.05227395   -0.7423603  -0.53745021
#> 12        4           2 -0.60215382       0.06820800   -0.7358390  -0.46846859
#> 10        4 (intercept)  0.47105039       0.12368367    0.2286349   0.71346592
#> 6         2           2 -0.13838624       0.03988563   -0.2165606  -0.06021185
#> 3         1           2 -0.04985085       0.04445023   -0.1369717   0.03727000
#>    Robust Wald p Robust Score p
#> 11  6.648178e-27             NA
#> 2   5.000849e-33             NA
#> 8   4.784250e-24             NA
#> 4   3.654085e-35             NA
#> 7   2.389181e-29             NA
#> 1   6.064411e-21             NA
#> 5   2.030015e-11             NA
#> 9   1.868292e-34             NA
#> 12  1.063748e-18             NA
#> 10  1.398072e-04             NA
#> 6   5.212790e-04             NA
#> 3   2.620759e-01             NA
weak_test$test_stat
#> [1] 10.02201
weak_test$p
#> [1] 0.006664209
```

We can see that the strong test has a larger robust Rao test statistic
and smaller p-value, compared to the weak test. This makes sense, as it
tests a stronger hypothesis that is more likely to be false (since it
contains all possible weak hypotheses).
