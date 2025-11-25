# Compute D matrix contribution for a given cluster to robust score statistic for glm robust score tests.

Compute D matrix contribution for a given cluster to robust score
statistic for glm robust score tests.

## Usage

``` r
D_matrix_contribution(indices, model_fits, yy, xx, family, link)
```

## Arguments

- indices:

  Indices of observations in cluster.

- model_fits:

  The fitted glm under the null hypothesis.

- yy:

  Vector of observed responses.

- xx:

  Design matrix for model.

- family:

  The model family for the fitted glm.

- link:

  The link function utilized in the fitted glm.
