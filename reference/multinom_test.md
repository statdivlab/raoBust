# Robust score (Rao) tests for multinomial regression.

Robust score (Rao) tests for multinomial regression.

## Usage

``` r
multinom_test(
  X = NULL,
  Y,
  formula = NULL,
  data = NULL,
  strong = FALSE,
  j = NULL,
  all_score = FALSE,
  penalty = FALSE,
  pseudo_inv = FALSE
)
```

## Arguments

- X:

  A \\n x p\\ design matrix of covariates.

- Y:

  A \\n x J\\ data matrix of outcomes.

- formula:

  a one-sided formula specifying the form of the mean model to be fit
  (use with `data` argument if `X` is not included)

- data:

  a dataframe with \\n\\ rows containing variables given in `formula`
  (use with `formula` argument if `X` is not included)

- strong:

  If FALSE, this function will compute the robust score statistic to
  test the weak null that for one specific \\j\\, \\\beta_j = 0\\ for
  the length \\p\\ vector \\\beta_j\\. If TRUE, this function instead
  computes the robust score statistic to test the strong null that
  \\\beta_1 = \beta_2 = \dots = \beta\_{J-1} = 0\\ for all length \\p\\
  vectors \\\beta_j\\, \\j\in\\1,\ldots,J-1\\\\. Default is FALSE.

- j:

  If `strong` is FALSE, this argument must be supplied. This gives the
  category \\j\\ in the weak null hypothesis that \\\beta_j = 0\\.

- all_score:

  If TRUE, score tests for each individual covariate and category pair
  (i.e. null that \\\beta\_{jk} = 0\\ for each category \\j = 1, \dots,
  J-1\\ and each covariate \\k = 1, \dots, p\\ pair) will be run and
  reported in output coefficient table. Default is FALSE.

- penalty:

  If TRUE will apply a Firth penalty to estimation under the alternative
  and under the null. Defaults to FALSE (ask Amy her preference)

- pseudo_inv:

  Use pseudo inverse for inverted portion of the robust score test to
  avoid issues with nearly singular matrices.

## Value

The robust score test statistic for the specified hypothesis test. A
list including the test statistic, p-value, estimated parameters under
the null hypothesis, and estimated parameters under the alternative
hypothesis.

## Author

Shirley Mathur
