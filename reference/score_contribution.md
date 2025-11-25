# Compute score contribution to robust score statistic for glm robust score tests.

Compute score contribution to robust score statistic for glm robust
score tests.

## Usage

``` r
score_contribution(i, model_fits, yy, xx, m = 1, family, link)
```

## Arguments

- i:

  Index of observation.

- model_fits:

  The fitted glm under the null hypothesis.

- yy:

  Vector of observed responses.

- xx:

  Design matrix for model.

- m:

  Number of parameters fixed under null hypothesis.

- family:

  The model family for the fitted glm.

- link:

  The link function utilized in the fitted glm.
