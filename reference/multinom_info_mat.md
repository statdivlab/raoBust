# Compute Information matrix for model parameters.

Compute Information matrix for model parameters.

## Usage

``` r
multinom_info_mat(X, Y, probs)
```

## Arguments

- X:

  This should be the \\n x p\\ design matrix of covariates.

- Y:

  This should be the \\n x J\\ data matrix of outcomes.

- probs:

  This should be the \\n \times J\\ matrix of estimated probabilities
  for each sample.

## Value

The Fisher information matrix for model parameters.

## Author

Shirley Mathur
