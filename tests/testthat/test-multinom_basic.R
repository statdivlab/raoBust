nsim <- 200
jj <- 5

test_that("beta vector to matrix function works as expected", {
  set.seed(240506)
  
  beta_true <- matrix(rpois(9,1), ncol = 3)
  null_j <- 2
  p <- nrow(beta_true) - 1
  J <- ncol(beta_true) + 1
  beta_values <- c(as.vector(beta_true[ ,-null_j]), beta_true[1,null_j])
  
  beta_copy <- beta_vector_to_matrix(beta_values, p, J, null_j)
  
  expect_true(beta_true[1,2] == beta_true[1,2])
  expect_true(beta_true[3,1] == beta_true[3,1])
  
})

test_that("multinomial test statistics are all positive", {

  nn <- 20
  set.seed(240506)
  ts_weak <- rep(NA, nsim)
  ts_strong <- rep(NA, nsim)

  df <- simulate_null_data_mult(nn = nn, strong=TRUE, sd_beta1s=2, sd_beta0s=1, jj_null=2)
  expect_silent({ result_strong <- get_multinom_score(df$X, df$Y, strong=TRUE) })
  expect_silent({ result_weak <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2) })

  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    result_weak <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
    result_strong <- get_multinom_score(df$X, df$Y, strong=TRUE)

    ts_weak[i] <- result_weak$test_stat
    ts_strong[i] <- result_strong$test_stat
  }
  
  expect_true(all(ts_weak > 0))
  expect_true(all(ts_strong > 0))

})

