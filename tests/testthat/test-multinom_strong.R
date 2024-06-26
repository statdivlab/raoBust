nsim <- 100
jj <- 5


test_that("strong multinomial test controls T1E under DGP for small n", {

  nn <- 20
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=TRUE)
    result <- multinom_test(df$X, df$Y, strong=TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))
  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10)) ## should be 1, 5, 10% to be nominal
  # for n = 500, 1000 sims: (3.6  6.4 10.5)%. High. Investigate.
  # for n = 10, 1000 sims: (3.4 6.0 7.1)%. Ok

})


test_that("strong multinomial test controls T1E under DGP for large n", {

  nn <- 500
  set.seed(240507)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=TRUE)
    result <- multinom_test(df$X, df$Y, strong=TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim)) ## empirical = 15%, upper cutpoint = 15.8%
  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))

})


test_that("strong multinomial test has increasing power under DGP with increasing n", {
  
  nn <- 15
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 100
  set.seed(240506)
  ts_big_n <- rep(NA, nsim)
  ps_big_n <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts_big_n[i] <- result$test_stat
    ps_big_n[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_n < 0.05))
  
})

test_that("strong multinomial test has increasing power under DGP with increasing magnitude", {
  
  nn <- 30
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 0.05, sd_beta1s = 0)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 30
  set.seed(240506)
  ts_big_b <- rep(NA, nsim)
  ps_big_b <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 2, sd_beta1s = 0)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts_big_b[i] <- result$test_stat
    ps_big_b[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_b < 0.05))
  
})

test_that("strong multinomial test controls T1E for overdispersed data for large n", {

  nn <- 500
  set.seed(2405071)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=TRUE, overdispersion=1)
    result <- multinom_test(df$X, df$Y, strong=TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("strong multinomial test controls T1E for overdispersed data for small n", {

  nn <- 20
  set.seed(2405072)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong = TRUE, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("strong multinomial test has increasing power for overdispersed data with increasing n", {
  
  nn <- 30
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 100
  set.seed(240506)
  ts_big_n <- rep(NA, nsim)
  ps_big_n <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts_big_n[i] <- result$test_stat
    ps_big_n[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_n < 0.05))
  
})

test_that("strong multinomial test has increasing power for overdispersed data with increasing magnitude", {
  
  nn <- 50
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 0.5, overdispersion = 1, sd_beta1s = 0)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 50
  set.seed(240506)
  ts_big_b <- rep(NA, nsim)
  ps_big_b <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 2, overdispersion = 1, sd_beta1s = 0)
    result <- multinom_test(df$X, df$Y, strong = TRUE)
    ts_big_b[i] <- result$test_stat
    ps_big_b[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_b < 0.05))
  
})
