nsim <- 200
jj <- 5

test_that("multinomial test statistics are all positive", {

  nn <- 20
  set.seed(240506)
  ts_weak <- rep(NA, nsim)
  ts_strong <- rep(NA, nsim)

  df <- simulate_null_data_mult(nn = nn, strong=TRUE, sd_beta1s=2, sd_beta0s=1, jj_null=2)
  expect_silent({ result_strong <- get_multinom_score(df$X, df$Y, strong=TRUE) })
  expect_silent({ result_weak <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2) })

  for (i in 1:79) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    result_weak <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
    result_strong <- get_multinom_score(df$X, df$Y, strong=TRUE)

    ts_weak[i] <- result_strong$test_stat
    ts_strong[i] <- result_weak$test_stat
  }

  expect_true(all(ts_weak > 0))
  expect_true(all(ts_strong > 0))

})

