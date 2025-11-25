# Run robust wald test for null hypothesis of form A x beta = c

Run robust wald test for null hypothesis of form A x beta = c

## Usage

``` r
lincom(test_object, A, c)
```

## Arguments

- test_object:

  Object of type `raoFit`.

- A:

  Matrix specifying linear combination of parameters for user-specified
  hypotheses.

- c:

  Vector that has user-specified null hypothesis value for
  user-specified linear combinations of parameters of interest.

## Value

Table with relevant quantities of interest for hypothesis tests.

## Author

Shirley Mathur

## Examples

``` r
#set true value of beta for DGP
beta0s <- rnorm(n = 4)
beta1s <- rnorm(n = 4)
beta_true <- rbind(beta0s, beta1s)
beta_true[,2] <- 0 
beta_true[2,1] <- 0 

#generate sample data
sample_data <- simulate_data_mult(30, Beta = beta_true)

#run weak multinom test
sample_multinom_test <- multinom_test(X = sample_data$X,
                                      Y = sample_data$Y,
                                      j = 2)

#test hypothesis that the coefficient for covariate is same for j = 1, j = 2
#so, we test hypothesis that beta_{k=1,j=1} - beta_{k=1, j=2} = 0

#first, set up A matrix
my_A <- set_up_lin_com(J = 5, p = 1, n_hypotheses = 1) 
my_A[1,"k_1_j_1"] <- 1
my_A[1, "k_1_j_2"] <- -1

#then, set c
my_c <- 0

#run linear combination test
sample_lincom_test <- lincom(sample_multinom_test,
                             A = my_A,
                             c = my_c)

#print result
sample_lincom_test
#>            Est   Std..Err Lower.CI..95.. Upper.CI..95..      pval
#> 1 -0.001795701 0.02477761    -0.05035892     0.04676752 0.9422257
```
