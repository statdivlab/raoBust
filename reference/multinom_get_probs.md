# Compute multinomial probabilities for a given value of model parameters.

Compute multinomial probabilities for a given value of model parameters.

## Usage

``` r
multinom_get_probs(X, Y, beta)
```

## Arguments

- X:

  This should be the \\n \times p\\ design matrix of covariates.

- Y:

  This should be the \\n \times J\\ data matrix of outcomes.

- beta:

  This should be the \\(p+1) \times (J-1) \beta\\ matrix.

## Value

The multinomial probabilities for a given value of \\\beta\\.

## Author

Shirley Mathur
