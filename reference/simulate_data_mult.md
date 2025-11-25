# Simulate multinomial data under the null (strong or weak) or alternative

By default, \\p = 2\\ with a single Uniform(0,1) continuous covariate.

## Usage

``` r
simulate_data_mult(
  nn,
  null = TRUE,
  strong = FALSE,
  alt_magnitude = 1,
  jj = 5,
  ms = 10000,
  jj_null = NULL,
  Beta = NULL,
  sd_beta0s = NULL,
  sd_beta1s = NULL,
  overdispersion = 0,
  covariate = NULL
)
```

## Arguments

- nn:

  Number of observations

- null:

  If TRUE will simulate under the null, if FALSE will simulate under the
  alternative

- strong:

  If null is TRUE, simulate under the strong null? Defaults to FALSE
  (simulate under the weak null)

- alt_magnitude:

  The mean of each parameter in the beta1 vector if null = FALSE.
  Defaults to \\1\\.

- jj:

  Number of taxa

- ms:

  Number of counts per sample

- jj_null:

  For the weak null, which taxon should be null

- Beta:

  User-specified value of the true beta (if you wish to draw from a
  fixed beta rather than generate beta0's and beta1's).

- sd_beta0s:

  The beta0's are drawn from a normal distribution with mean zero. This
  is the standard deviation of that distribution.

- sd_beta1s:

  The beta1's are drawn from a normal distribution with mean zero under
  the null or non-zero under the alternative. This is the standard
  deviation of that distribution.

- overdispersion:

  An additional normal random variable can be added to the link function
  to add dispersion above a multinomial distribution. This is the
  standard aviation for that normal variable. Useful for confirming
  error rate control under model misspecification.

- covariate:

  An optional covariate vector, if not provided the covariate will be a
  sequence of increasing values from 0 to 1.

## Author

Amy Willis
