test_that("there are large estimates without Firth penalty under separation", {
  
  df <- simulate_data_mult(nn = 100, strong = TRUE, covariate = rep(0:1, each = 50))
  df$Y[, 1] <- 0
  df$Y[1:50, 1] <- 1000
  
  result <- multinom_test(df$X, df$Y, strong = TRUE)
  result_penalty <- multinom_test(df$X, df$Y, strong = TRUE, penalty = TRUE)
  # expect that B_1 for separated category will be largest parameter 
  expect_true(which.max(abs(result$mle1)) == 2)
  # expect that when adding the Firth penalty the B_1 will now have a smaller magnitude
  expect_true(abs(result$mle1[2, 1]) > abs(result_penalty$mle1[2, 1]))
  
})

test_that("Firth penalized estimates are close to those from radEmu", {
  
  skip("Skipping to avoid including radEmu as a Suggests, can check this manually")
  
  df <- simulate_data_mult(nn = 100, null = FALSE, covariate = rep(0:1, each = 5))
  result <- multinom_test(df$X, df$Y, strong = FALSE, j = 2, penalty = TRUE)
  
  lcc <- function(x) { return(x[length(x)]) }
  lcc_grad <- function(x) {
    vec <- rep(0, length = length(x)) 
    vec[length(x)] <- 1
    return(vec)
  }
  result_radEmu <- emuFit(Y = df$Y, X = cbind(1, df$X), run_score_tests = TRUE,
                          compute_cis = FALSE, constraint_fn = lcc,
                          constraint_grad_fn = lcc_grad, constraint_param = NA,
                          test_kj = data.frame(k = 2, j = 2), return_nullB = TRUE)
  # check on MLE under alternative
  expect_true(all.equal(as.vector(result$mle1), 
                        as.vector(result_radEmu$B[, -5]), 
                        tolerance = 0.01*max(abs(result$mle1))))
  # check on MLE under null
  expect_true(all.equal(as.vector(result$mle0), 
                        as.vector(result_radEmu$null_B[[1]][, -5]), 
                        tolerance = 0.01*max(abs(result$mle0))))
  
})


test_that("Firth penalized estimates are close to those from radEmu", {
  
  skip("Skipping to avoid including radEmu as a Suggests, can check this manually")
  
  df <- simulate_data_mult(nn = 100, null = FALSE, covariate = rep(0:1, each = 5))
  result <- multinom_test(df$X, df$Y, strong = TRUE, penalty = TRUE)
  
  lcc <- function(x) { return(x[length(x)]) }
  lcc_grad <- function(x) {
    vec <- rep(0, length = length(x)) 
    vec[length(x)] <- 1
    return(vec)
  }
  result_radEmu <- emuFit(Y = df$Y, X = cbind(1, df$X), run_score_tests = FALSE,
                          compute_cis = FALSE, constraint_fn = lcc,
                          constraint_grad_fn = lcc_grad, constraint_param = NA)
  expect_true(all.equal(as.vector(result$mle1), 
                        as.vector(result_radEmu$B[, -5]), tolerance = 0.01))
  
})

test_that("Firth penalized estimates are close to unpenalized estimates for large n", {
  df <- simulate_data_mult(nn = 100, null = FALSE, covariate = rep(0:1, each = 5))
  result <- multinom_test(df$X, df$Y, strong = TRUE, penalty = TRUE)
  result_np <- multinom_test(df$X, df$Y, strong = TRUE)
  # check for estimation under alternative
  expect_true(all.equal(as.vector(result$mle1), 
                        as.vector(result_np$mle1), tolerance = 0.01))
  # check for estimation under null
  expect_true(all.equal(as.vector(result$mle0), 
                        as.vector(result_np$mle0), tolerance = 0.01))
})
