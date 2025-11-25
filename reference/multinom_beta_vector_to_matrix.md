# Create \\\beta\\ matrix from vector of \\\beta\\ for \\\beta_k : k \neq j\\

Create \\\beta\\ matrix from vector of \\\beta\\ for \\\beta_k : k \neq
j\\

## Usage

``` r
multinom_beta_vector_to_matrix(values, p, J, null_j, beta_j_null = NULL)
```

## Arguments

- values:

  A vector containing the values for all \\\beta_k\\, with \\k \neq j\\,
  as well as all \\\beta\_{k0}, for k = 1, \dots, J\\. In particular,
  this vector should be so that the first \\(J-2)(p+1)\\ entries are
  \\\beta\_{10}, \beta\_{1}^{\top}, \beta\_{20}, \beta\_{2}^{\top},
  \dots, \beta\_{(j-1)0}, \beta\_{j-1}^{\top}, \beta\_{(j+1)0},
  \beta\_{j+1}^{\top}, \dots, \beta\_{(J-1)0}, \beta\_{J-1}^{\top},
  \beta\_{j0}\\. Then, the \\(J-2)(p+1) + 1\\ entry should be
  \\\beta\_{j0}\\.

- p:

  This should be the number of covariates.

- J:

  This should be the number of categories.

- null_j:

  This specifies for which category you have set \\\beta_j = 0\\.

- beta_j_null:

  This is the null hypothesized value for \\\beta_j\\, which is by
  default set to be \\\beta_j = 0\\.

## Value

The full \\(p+1) \times (J-1)\\ matrix of \\\beta\\

## Author

Shirley Mathur
