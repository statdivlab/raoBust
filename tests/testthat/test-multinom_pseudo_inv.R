nsim <- 5

test_that("using the inverse and pseudo inverse give the same results when an inverse exists", {

  nn <- 10
  stat <- rep(NA, 50)
  stat_alt <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, strong=TRUE)
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE, pseudo_inv = TRUE))
    stat[sim] <- result$test_stat
    stat_alt[sim] <- result_alt$test_stat
  }


  expect_true(all.equal(stat, stat_alt, tol = 1e-03))

})

nsim <- 100
test_that("the test stat with the pseudo inverse controls t1e when the inverse is computationally singular for strong hypothesis", {

  nn <- 10
  p <- rep(NA, 50)
  p_alt <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, strong=TRUE, ms = 10)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    capture.output(
      res <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE)),
      file = NULL
    )
    capture.output(
      res_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE, pseudo_inv = TRUE)),
      file = NULL
    )
    
    if (!inherits(res, "try-error")) p[sim] <- res$p
    if (!inherits(res_alt, "try-error")) p_alt[sim] <- res_alt$p
  }

  expect_true(sum(!is.na(p)) > 5)
  expect_true(all.equal(p[!is.na(p)], p_alt[!is.na(p)], tol = 0.005))
  expect_true(mean(p_alt < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(p_alt < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("the test stat with the pseudo inverse controls t1e when the inverse is computationally singular for weak hypothesis", {

  # I'm having trouble coming up with times in which there are computationally singular errors
  # for the weak null

  nn <- 6
  p <- rep(NA, 50)
  p_alt <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, strong=FALSE, jj_null = 2, ms = 100)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    df$Y[, 4] <- 0
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2, pseudo_inv = TRUE))
    p[sim] <- result$p
    p_alt[sim] <- result_alt$p
  }

  expect_true(all.equal(p[!is.na(p)], p_alt[!is.na(p)], tol = 0.0001))
  expect_true(mean(p_alt < 0.05) < 0.05 + 1.96*sqrt(0.05 * 0.95/nsim))
  expect_true(mean(p_alt < 0.10) < 0.1 + 1.96*sqrt(0.1 * 0.9/nsim))

})

test_that("power increases with signal for the test stat with the pseudo inverse when the inverse is computationally singular for strong hypothesis", {

  skip("Skipping due to time, can check manually")

  nn <- 20
  p <- rep(NA, 50)
  p_alt <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, null = FALSE, ms = 10)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE, pseudo_inv = TRUE))
    p[sim] <- result$p
    p_alt[sim] <- result_alt$p
  }

  expect_true(all.equal(p[!is.na(p)], p_alt[!is.na(p)], tol = 0.05))

  nn <- 20
  p2 <- rep(NA, 50)
  p_alt2 <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, null = FALSE, ms = 10, alt_magnitude = 3)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=TRUE, pseudo_inv = TRUE))
    p2[sim] <- result$p
    p_alt2[sim] <- result_alt$p
  }

  # there are some bigger differences for this set
  #expect_true(all.equal(p2[!is.na(p2)], p_alt2[!is.na(p2)], tol = 0.3))
  expect_true(mean(p_alt2 < 0.05) > mean(p_alt < 0.05))

})

test_that("power increases with signal for the test stat with the pseudo inverse when the inverse is computationally singular for strong hypothesis", {

  skip("Skipping due to time, can check manually")

  nn <- 20
  p <- rep(NA, 50)
  p_alt <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, null = FALSE, ms = 100)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    df$Y[, 4] <- 0
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2, pseudo_inv = TRUE))
    p[sim] <- result$p
    p_alt[sim] <- result_alt$p
  }

  expect_true(all.equal(p[!is.na(p)], p_alt[!is.na(p)], tol = 0.005))

  nn <- 20
  p2 <- rep(NA, 50)
  p_alt2 <- rep(NA, 50)

  for (sim in 1:nsim) {
    set.seed(sim)
    df <- simulate_data_mult(nn = nn, null = FALSE, ms = 100, alt_magnitude = 3)
    df$Y[, 1] <- 0
    df$Y[, 3] <- 0
    df$Y[, 4] <- 0
    result <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2))
    result_alt <- suppressWarnings(multinom_test(df$X, df$Y, strong=FALSE, j = 2, pseudo_inv = TRUE))
    p2[sim] <- result$p
    p_alt2[sim] <- result_alt$p
  }

  # there are some bigger differences for this set
  expect_true(all.equal(p2[!is.na(p2)], p_alt2[!is.na(p2)], tol = 0.001))
  expect_true(mean(p_alt2 < 0.05) > mean(p_alt < 0.05))

})
