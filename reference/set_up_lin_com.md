# Create matrix to be populated with coefficients for user-specified number of hypotheses.

Create matrix to be populated with coefficients for user-specified
number of hypotheses.

## Usage

``` r
set_up_lin_com(J, p, n_hypotheses)
```

## Arguments

- J:

  Number of categories.

- p:

  Number of coefficients (excluding intercept).

- n_hypotheses:

  Number of hypotheses to test.

## Value

An 'A' matrix with `n_hypotheses` rows and `(J-1) x (p+1)` columns with
the various coefficient and category combinations as column names. The
matrix is filled with all 0's, but can be modified by the user to
reflect the linear combinations they wish to test.

## Author

Shirley Mathur
