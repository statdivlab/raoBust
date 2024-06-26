test_that("there are large estimates without Firth penalty under separation", {
  
  df <- simulate_data_mult(nn = 10, strong = TRUE)
  df$X <- matrix(rep(0:1, each = 5), nrow = 10, ncol = 1)
  df$Y[, 1] <- 0
  df$Y[1, 1] <- 1000
  
  result <- multinom_test(df$X, df$Y, strong = TRUE)
  
})
