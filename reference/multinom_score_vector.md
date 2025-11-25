# Compute score for model parameters.

Compute score for model parameters.

## Usage

``` r
multinom_score_vector(X, Y, probs)
```

## Arguments

- X:

  This should be the \\n \times p\\ design matrix of covariates.

- Y:

  This should be the \\n \times J\\ data matrix of outcomes.

- probs:

  This should be the \\n \times J\\ matrix of estimated probabilities
  for each sample.

## Value

The score vector for the model parameters.

## Author

Shirley Mathur
