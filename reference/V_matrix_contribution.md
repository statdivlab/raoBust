# Compute V matrix contribution to robust score statistic for glm robust score tests.

Compute V matrix contribution to robust score statistic for glm robust
score tests.

## Usage

``` r
V_matrix_contribution(
  indices,
  model_fits,
  yy,
  xx,
  pp0 = 1,
  corr_mat,
  family,
  link
)
```

## Arguments

- indices:

  Indices of observations for cluster.

- model_fits:

  The fitted glm under the null hypothesis for the given indices.

- yy:

  Vector of observed responses.

- xx:

  Design matrix for model.

- pp0:

  Number of fixed parameters under null hypothesis.

- corr_mat:

  Working correlation matrix for model fit.

- family:

  The model family for the fitted glm.

- link:

  The link function utilized in the fitted glm.
