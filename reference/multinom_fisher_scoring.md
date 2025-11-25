# Optimization under null or alternative for multinomial model via Fisher scoring.

Optimization under null or alternative for multinomial model via Fisher
scoring.

## Usage

``` r
multinom_fisher_scoring(
  beta,
  X,
  Y,
  null = TRUE,
  strong = FALSE,
  null_j = NULL,
  j_ind = NULL,
  k_ind = NULL,
  tol = 1e-05,
  stepSize = 0.5,
  arm_c = 0.5,
  maxit = 250,
  pseudo_inv = FALSE
)
```

## Arguments

- beta:

  The initial values provided for the \\\beta\\ parameters.

- X:

  The \\n x p\\ design matrix of covariates.

- Y:

  The \\n x J\\ data matrix of outcomes.

- null:

  If TRUE, optimizes under the null, if FALSE, optimizes under the
  alternative. Defaults to TRUE.

- strong:

  If FALSE, this function will compute the robust score statistic to
  test the weak null that for one specific \\j\\, \\\beta_j = 0\\ for
  the length \\p\\ vector \\\beta_j\\. If TRUE, this function instead
  computes the robust score statistic to test the strong null that
  \\\beta_1 = \beta_2 = \dots = \beta\_{J-1} = 0\\ for all length \\p\\
  vectors \\\beta_j\\, \\j\in\\1,\ldots,J-1\\\\. Default is FALSE.

- null_j:

  If `strong` is FALSE, this argument must be supplied. This gives the
  category \\j\\ in the weak null hypothesis that \\\beta_j = 0\\.
  Default is NULL.

- j_ind:

  If `strong` is FALSE and `null_j` is NULL, this argument must be
  supplied. This gives the category index of the individual covariate
  that is tested in the weak null hypothesis that \\\beta\_{kj} = 0\\.

- k_ind:

  If `strong` is FALSE and `null_j` is NULL, this argument must be
  supplied. This gives the covariate index of the individual covariate
  that is tested in the weak null hypothesis that \\\beta\_{kj} = 0\\.

- tol:

  The tolerance used to determine how much better update function value
  must be prior to stopping algorithm.

- stepSize:

  The size of the step to take during the parameter update step.

- arm_c:

  Control parameter for checking Armijo condition.

- maxit:

  Maximum number of iterations for Fisher scoring. Defaults to 250.

- pseudo_inv:

  Use the pseudo-inverse of the Fisher information matrix for the update
  (in case the inverse in computationally singular)

## Value

The optimal beta values under the null or alternative model.

## Author

Shirley Mathur
