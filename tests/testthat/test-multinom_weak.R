nsim <- 200
jj <- 5



test_that("weak multinomial test controls T1E under DGP for small n", {

  nn <- 20
  set.seed(240506)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    ts[i] <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
  }

  ps <- pchisq(ts, df = 1, lower.tail=FALSE)
  hist(ps)

  100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})



test_that("weak multinomial test controls T1E under DGP for large n", {

  nn <- 1000
  set.seed(240506)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    ts[i] <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
  }

  ps <- pchisq(ts, df = 1, lower.tail=FALSE)
  hist(ps)

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))


})




test_that("weak multinomial test controls T1E for overdispersed data for small n", {

  nn <- 20
  set.seed(2405063)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2, overdispersion=2)
    ts[i] <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
  }

  ps <- pchisq(ts, df = 1, lower.tail=FALSE)
  hist(ps)

  100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})



test_that("weak multinomial test controls T1E for overdispersed data for large n", {

  nn <- 200
  set.seed(2405064)
  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2, overdispersion=2)
    ts[i] <- get_multinom_score(df$X, df$Y, strong=FALSE, j=2)
  }

  ps <- pchisq(ts, df = 1, lower.tail=FALSE)
  hist(ps)

  100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))


})


