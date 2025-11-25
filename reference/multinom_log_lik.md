# Negative log-likelihood for multinomial data under the alternative

Negative log-likelihood for multinomial data under the alternative

## Usage

``` r
multinom_log_lik(beta_as_vector, Y, X)
```

## Arguments

- beta_as_vector:

  A vector containing the values for all \\\beta_k\\ for \\k = 1, \dots,
  J-1\\ and \\\beta\_{k0}, for k = 1, \dots, J\\. In particular, this
  vector should be so that the \\(J-1)(p+1)\\ entries are \\\beta\_{10},
  \beta\_{1}^{\top}, \beta\_{20}, \beta\_{2}^{\top}, \dots,
  \beta\_{(J-1)0}, \beta\_{J-1}^{\top}, \beta\_{j0}\\.

- Y:

  This should be the \\n x J\\ data matrix of outcomes.

- X:

  This should be the \\n x p\\ design matrix of covariates.

## Value

The value of the log likelihood for the input \\\beta\\

## Author

Shirley Mathur
