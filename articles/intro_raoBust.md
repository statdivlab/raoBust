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
link), binomial GLMs (logit link), linear models (identity link), and
multinomial GLMs (log link).

## Model fitting

### GLMs

We will demonstrate the function
[`glm_test()`](../reference/glm_test.md) in `raoBust` using the
`light_presence` dataset.

``` r
data("light_presence")
head(light_presence)
#> # A tibble: 6 × 5
#>   sample   cog     presence clade light
#>   <chr>    <chr>      <dbl> <chr> <chr>
#> 1 AS9601   COG4251        1 HL_II HL   
#> 2 CCMP1375 COG4251        0 LL_II LL   
#> 3 EQPAC1   COG4251        1 HL_I  HL   
#> 4 GP2      COG4251        1 HL_II HL   
#> 5 LG       COG4251        0 LL_II LL   
#> 6 MED4     COG4251        1 HL_I  HL
```

This dataset provides presence/absence data for COG 4251 with a set of
*Prochlorococcus* genomes. COG 4251 has the function Bacteriophytochrome
(light-regulated signal transduction histidine kinase). We would like to
test whether the presence or absence of this COG is associated with a
genome being from a high light or low light condition. We will use a
binomial GLM model to investigate this.

``` r
glm_test(formula = presence ~ light, data = light_presence, family = binomial(link = "logit"))
#> 
#> Call:
#> glm(formula = presence ~ light, data = light_presence, family = binomial(link = "logit"))
#> 
#> 
#> Coefficient estimates:
#>              Estimate Non-robust Std Error Robust Std Error Lower 95% CI
#> (Intercept)  20.56607             3964.631     7.222006e-09     20.56607
#> lightLL     -20.74839             3964.631     6.260623e-01    -21.97545
#>             Upper 95% CI Non-robust Wald p Robust Wald p Robust Score p
#> (Intercept)     20.56607         0.9958611             0   5.465931e-06
#> lightLL        -19.52133         0.9958244             0   3.327047e-03
```

From the output, we can see that we have an estimated log odds of COG
4251 being present in a low light compared to a high light condition of
$- 20.75$, with a robust Rao test p-value of $0.003$. This means that
COG 4251 is present far less often in low light compared to high light.
This makes sense if we look at the presence/absence breakdown for this
COG.

``` r
table(light_presence[, c("presence", "light")])
#>         light
#> presence HL LL
#>        0  0  6
#>        1 20  5
```

We can see that this COG is present in all samples from the high light
condition, and only in 5 out of 11 samples from the low light condition.

We can compare the results from [`glm_test()`](../reference/glm_test.md)
to the results we would get if we ran a typical
[`glm()`](https://rdrr.io/r/stats/glm.html).

``` r
summary(glm(formula = presence ~ light, data = light_presence, family = binomial(link = "logit")))
#> 
#> Call:
#> glm(formula = presence ~ light, family = binomial(link = "logit"), 
#>     data = light_presence)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)    20.57    3964.63   0.005    0.996
#> lightLL       -20.75    3964.63  -0.005    0.996
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 30.462  on 30  degrees of freedom
#> Residual deviance: 15.158  on 29  degrees of freedom
#> AIC: 19.158
#> 
#> Number of Fisher Scoring iterations: 19
```

Here we see the same estimates, but `glm` does not give us robust
standard errors or robust Wald and Rao tests.

### GEEs, accounting for clusters

Next, we will use the function [`gee_test()`](../reference/gee_test.md)
to account for clustered observations. We will apply this to the `ddPCR`
dataset, which gives digital droplet PCR measurements from participants
in a study of people living with cystic fibrosis.

``` r
data("ddPCR")
#> Warning in data("ddPCR"): data set 'ddPCR' not found
head(ddPCR)
#>   Subject Treatment SpecimenType   Value
#> 1       1         0       Sputum 2457600
#> 2       1         0       Sputum 2081314
#> 3       1         0       Sputum 2554822
#> 4       1         0       Sputum 3122037
#> 5       1         1       Sputum  855946
#> 6       2         0       Sputum  647279
```

We will fit a Poisson GEE model to estimate and test the expected
fold-difference in the `ddPCR` measurement associated with whether or
not the participant is receiving the treatment, accounting for the
specimen type (specimens are collected either from sputum or saliva). We
choose Poisson regression because we want to estimate a fold-difference
in means, and because we are using robust tests we will not worry about
whether we actually expect `ddPCR` measurements to follow a Poisson
distribution. We use a GEE instead of a GLM because we have multiple
observations from the same samples.

``` r
gee_test(formula = Value ~ Treatment + SpecimenType, data = ddPCR, family = poisson(link = "log"), id = Subject)
#> 
#> Call:
#> gee_test(formula = Value ~ Treatment + SpecimenType, data = ddPCR, 
#>     family = poisson(link = "log"), id = Subject)
#> 
#> 
#> Coefficient estimates:
#>                       Estimate Robust Std Error Lower 95% CI Upper 95% CI
#> (Intercept)        13.78127232        0.2355520  13.31959889   14.2429458
#> Treatment          -0.30102500        0.3629500  -1.01239389    0.4103439
#> SpecimenTypeSputum  0.73510333        0.3532041   0.04283611    1.4273706
#> correlation:alpha   0.04471394               NA           NA           NA
#>                    Robust Wald p Robust Score p
#> (Intercept)            0.0000000     0.01928004
#> Treatment              0.4068870     0.44146433
#> SpecimenTypeSputum     0.0374117     0.09694808
#> correlation:alpha             NA             NA
```

We estimate an expected fold-difference of $\exp( - 0.32) = 0.72$ in the
`ddPCR` measurement for participants who have received the treatment,
accounting for sample type, with a robust Rao p-value of $0.44$. We can
compare these results to what we would get if we ran a GLM, ignoring
correlation of samples from the same subjects.

``` r
glm_test(formula = Value ~ Treatment + SpecimenType, data = ddPCR, family = poisson(link = "log"))
#> 
#> Call:
#> glm(formula = Value ~ Treatment + SpecimenType, data = ddPCR, 
#>     family = poisson(link = "log"))
#> 
#> 
#> Coefficient estimates:
#>                      Estimate Non-robust Std Error Robust Std Error
#> (Intercept)        13.7883695         0.0001742993        0.2262519
#> Treatment          -0.3185865         0.0002952090        0.3225658
#> SpecimenTypeSputum  0.7335946         0.0002073782        0.2827819
#>                    Lower 95% CI Upper 95% CI Non-robust Wald p Robust Wald p
#> (Intercept)          13.3449239   14.2318151                 0   0.000000000
#> Treatment            -0.9508038    0.3136308                 0   0.323317364
#> SpecimenTypeSputum    0.1793523    1.2878369                 0   0.009480956
#>                    Robust Score p
#> (Intercept)          0.0002177321
#> Treatment            0.2943003843
#> SpecimenTypeSputum   0.0129525938
```

We can see that we get the same estimates as in the
[`gee_test()`](../reference/gee_test.md), but different standard errors
and p-values because we have not appropriately accounted for the
clustering.
