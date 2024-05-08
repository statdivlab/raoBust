nsim <- 100
jj <- 5


test_that("strong multinomial test controls T1E under DGP for small n", {

  nn <- 20
  set.seed(240506)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE)
    ts[i] <- get_multinom_score(df$X, df$Y, strong= TRUE)
  }
  ts[ts < 0] <- 0

  ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)

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
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE)
    ts[i] <- get_multinom_score(df$X, df$Y, strong= TRUE)
  }
  ts[ts < 0] <- 0

  ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)

  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim)) ## empirical = 15%, upper cutpoint = 15.8%
  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))

})


test_that("strong multinomial test controls T1E for overdispersed data for large n", {

  nn <- 500
  set.seed(2405071)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE, overdispersion=1)
    ts[i] <- get_multinom_score(df$X, df$Y, strong= TRUE)
  }
  ts[ts < 0] <- 0

  ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("strong multinomial test controls T1E for overdispersed data for small n", {

  nn <- 20
  set.seed(2405072)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=TRUE, overdispersion=1)
    ts[i] <- get_multinom_score(df$X, df$Y, strong= TRUE)
  }
  ts[ts < 0] <- 0

  ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})
