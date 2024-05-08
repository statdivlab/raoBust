nsim <- 200
jj <- 5

test_that("weak multinomial test statistics are all positive", {

  nn <- 20
  set.seed(240506)
  ts_weak <- rep(NA, nsim)
  ts_strong <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    ts_weak[i] <- get_multinom_score(df$X, df$Y, strong=TRUE)
    ts_strong[i] <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
  }

  expect_true(all(ts_weak > 0))
  expect_true(all(ts_strong > 0))

})

