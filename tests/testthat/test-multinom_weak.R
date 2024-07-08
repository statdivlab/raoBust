nsim <- 100
jj <- 5


test_that("weak multinomial test controls T1E under DGP for small n", {

  nn <- 20
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong = FALSE, sd_beta1s = 2, sd_beta0s = 0, jj_null = 2)
    result <- multinom_test(df$X, df$Y, strong = FALSE, j = 2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  #100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("weak multinomial test controls T1E under DGP for large n", {

  nn <- 1000
  set.seed(2405062)
  #set.seed(234738)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})


test_that("weak multinomial test has increasing power under DGP with increasing n", {
  
  nn <- 10
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 100
  set.seed(240506)
  ts_big_n <- rep(NA, nsim)
  ps_big_n <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts_big_n[i] <- result$test_stat
    ps_big_n[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_n < 0.05))
  
})

test_that("weak multinomial test has increasing power under DGP with increasing magnitude", {
  
  nn <- 50
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 0.5)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 50
  set.seed(240506)
  ts_big_b <- rep(NA, nsim)
  ps_big_b <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 2)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts_big_b[i] <- result$test_stat
    ps_big_b[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_b < 0.05))
  
})

test_that("weak multinomial test controls T1E for overdispersed data for small n", {

  nn <- 20
  set.seed(2405063)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2, overdispersion=2)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

 #  hist(ps)
  100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})



test_that("weak multinomial test controls T1E for overdispersed data for large n", {


  nn <- 200
  set.seed(2405064)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=0, jj_null=2, overdispersion=2)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  expect_true(mean(ps < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))


})

test_that("weak multinomial test has increasing power with overdispersed with increasing n", {
  
  nn <- 10
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 100
  set.seed(240506)
  ts_big_n <- rep(NA, nsim)
  ps_big_n <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts_big_n[i] <- result$test_stat
    ps_big_n[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_n < 0.05))
  
})

test_that("weak multinomial test has increasing power under DGP with increasing magnitude", {
  
  nn <- 50
  set.seed(240506)
  ts <- rep(NA, nsim)
  ps <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 0.5, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts[i] <- result$test_stat
    ps[i] <- result$p
  }
  
  nn <- 50
  set.seed(240506)
  ts_big_b <- rep(NA, nsim)
  ps_big_b <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, null = FALSE, alt_magnitude = 2, overdispersion = 1)
    result <- multinom_test(df$X, df$Y, strong=FALSE, j=2)
    ts_big_b[i] <- result$test_stat
    ps_big_b[i] <- result$p
  }
  
  expect_true(mean(ps < 0.05) < mean(ps_big_b < 0.05))
  
})




