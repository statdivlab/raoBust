nsim <- 200
jj <- 5

test_that("beta vector to matrix function works as expected", {
  set.seed(240506)
  
  beta_true <- matrix(rpois(9, 1), ncol = 3)
  null_j <- 2
  p <- nrow(beta_true) - 1
  J <- ncol(beta_true) + 1
  beta_values <- c(as.vector(beta_true[, -null_j]), beta_true[1, null_j])
  
  beta_copy <- beta_vector_to_matrix(beta_values, p, J, null_j)
  
  expect_true(beta_true[1, 2] == beta_copy[1, 2])
  expect_true(beta_true[3, 1] == beta_copy[3, 1])
  
})

test_that("simulate_data_mult function works as expected", {
  nn <- 20
  
  B <- simulate_data_mult(nn = nn, null = TRUE, jj_null = 2)$B
  expect_true(B[2, 2] == 0)
  
  B <- simulate_data_mult(nn = nn, null = TRUE, strong = TRUE, jj_null = 2)$B
  expect_true(all.equal(B[2, ], rep(0, jj - 1)))
  
  B <- simulate_data_mult(nn = nn, null = FALSE, jj_null = 2, alt_magnitude = 5)$B
  expect_true(B[2, 2] != 0)
  
})

test_that("errors and warnings work as expected", {
  nn <- 20
  df <- simulate_data_mult(nn = nn, strong = TRUE, sd_beta1s = 2, sd_beta0s = 1, jj_null = 2)
  expect_error(multinom_test(df$X, df$Y, strong = FALSE))
  
  df2 <- df
  df2$X <- df2$X[1:(nrow(df$X) - 1), , drop = FALSE]
  expect_error(multinom_test(df2$X, df2$Y, strong = FALSE, j = 2))
  
  df3 <- df
  rownames(df3$X) <- paste0("samp", 1:nrow(df3$X))
  rownames(df3$Y) <- paste0("obs", 1:nrow(df3$Y))
  expect_warning({res <- multinom_test(df3$X, df3$Y, strong = FALSE, j = 2)})
})

test_that("multinomial test statistics are all positive", {

  nn <- 20
  set.seed(240506)
  ts_weak <- rep(NA, nsim)
  ts_strong <- rep(NA, nsim)

  df <- simulate_data_mult(nn = nn, strong = TRUE, sd_beta1s = 2, sd_beta0s = 1, jj_null = 2)
  expect_silent({ result_strong <- multinom_test(df$X, df$Y, strong=TRUE) })
  expect_silent({ result_weak <- multinom_test(df$X, df$Y, strong=FALSE, j=2) })

  for (i in 1:nsim) {
    df <- simulate_data_mult(nn = nn, strong = TRUE, sd_beta1s = 2, sd_beta0s = 0, jj_null = 2)
    result_weak <- multinom_test(df$X, df$Y, strong = FALSE, j = 2)
    result_strong <- multinom_test(df$X, df$Y, strong = TRUE)

    ts_weak[i] <- result_weak$test_stat
    ts_strong[i] <- result_strong$test_stat
  }
  
  expect_true(all(ts_weak > 0))
  expect_true(all(ts_strong > 0))

})

test_that("estimation under the alternative is not strictly worse than nnet", {
  
  nsim <- 10
  nn <- 500
  
  neg_ll_df <- data.frame(raoBust = rep(NA, nsim),
                          nnet = rep(NA, nsim))
  
  for (sim in 1:nsim) {
    set.seed(sim)
    
    df <- simulate_data_mult(nn = nn, strong = TRUE, sd_beta1s = 2, sd_beta0s = 1, jj_null = 2)
    
    # estimate with raoBust
    result_strong5 <- multinom_test(df$X, df$Y, strong = TRUE)
    
    # estimate with nnet 
    result_nnet5 <- nnet::multinom(cbind(df$Y[, 5], df$Y[, 1:4]) ~ df$X)
    result_nnet5_summ <- summary(result_nnet5)
    
    # negative log likelihood 
    neg_ll_df[sim, 1] <- multinom_log_lik_alternative(as.vector(result_strong5$mle1), df$Y, df$X)
    neg_ll_df[sim, 2] <- multinom_log_lik_alternative(as.vector(t(result_nnet5_summ$coefficients)), df$Y, df$X)
   
    expect_true(all.equal(as.vector(result_strong5$mle1), as.vector(t(result_nnet5_summ$coefficients)), 
                          tolerance = 0.05)) 
  }
  
  # expect that raoBust won't be much worse than nnet
  expect_true(all.equal(neg_ll_df$raoBust, neg_ll_df$nnet, tolerance = 0.01*max(neg_ll_df$raoBust)))
  
})

test_that("estimation under the alternative is not strictly worse than nnet, data generated under alternative", {
  
  nsim <- 10
  nn <- 500
  
  neg_ll_df <- data.frame(raoBust = rep(NA, nsim),
                          nnet = rep(NA, nsim))
  
  for (sim in 1:nsim) {
    set.seed(sim)
    
    df <- simulate_data_mult(nn = nn, null = FALSE, sd_beta1s = 2, sd_beta0s = 1, jj_null = 2)
    
    # estimate with raoBust
    result_strong5 <- multinom_test(df$X, df$Y, strong = TRUE)
    
    # estimate with nnet 
    result_nnet5 <- nnet::multinom(cbind(df$Y[, 5], df$Y[, 1:4]) ~ df$X)
    result_nnet5_summ <- summary(result_nnet5)
    
    # negative log likelihood 
    neg_ll_df[sim, 1] <- multinom_log_lik_alternative(as.vector(result_strong5$mle1), df$Y, df$X)
    neg_ll_df[sim, 2] <- multinom_log_lik_alternative(as.vector(t(result_nnet5_summ$coefficients)), df$Y, df$X)
    
    expect_true(all.equal(as.vector(result_strong5$mle1), as.vector(t(result_nnet5_summ$coefficients)), 
                 tolerance = 0.05))
  }
  
  # expect that raoBust won't be much worse than nnet
  expect_true(all.equal(neg_ll_df$raoBust, neg_ll_df$nnet, tolerance = 0.01*max(neg_ll_df$raoBust)))
  
})

test_that("estimation under the strong null is not strictly worse than nnet and matches nnet pretty well", {
  
  nsim <- 10
  nn <- 500
  
  neg_ll_df <- data.frame(raoBust = rep(NA, nsim),
                          nnet = rep(NA, nsim))
  
  for (sim in 1:nsim) {
    set.seed(sim)
    
    df <- simulate_data_mult(nn = nn, strong = TRUE, sd_beta1s = 2, sd_beta0s = 1, jj_null = 2)
    
    # estimate with raoBust
    result_strong5 <- multinom_test(df$X, df$Y, strong = TRUE)
    
    # estimate with nnet 
    result_nnet5 <- nnet::multinom(cbind(df$Y[, 5], df$Y[, 1:4]) ~ 1)
    result_nnet5_summ <- summary(result_nnet5)
    
    # negative log likelihood 
    neg_ll_df[sim, 1] <- multinom_log_lik_alternative(as.vector(result_strong5$mle0), df$Y, df$X)
    neg_ll_df[sim, 2] <- multinom_log_lik_alternative(as.vector(t(result_nnet5_summ$coefficients)), df$Y, df$X)

    expect_equal(result_strong5$mle0[1, ], as.vector(result_nnet5_summ$coefficients), 
                 tolerance = 0.05)
  }
  
  summary(neg_ll_df$raoBust - neg_ll_df$nnet)
  
  # expect that nnet won't have a lower negative ll than raoBust every time
  expect_false(sum(neg_ll_df$raoBust - neg_ll_df$nnet < 0) == 0)
  
})

test_that("estimation under the strong null is not strictly worse than nnet and matches nnet pretty well, data generated under alternative", {
  
  nsim <- 10
  nn <- 500
  
  neg_ll_df <- data.frame(raoBust = rep(NA, nsim),
                          nnet = rep(NA, nsim))
  
  for (sim in 1:nsim) {
    set.seed(sim)
    
    df <- simulate_data_mult(nn = nn, null = FALSE, sd_beta1s=2, sd_beta0s=1, jj_null=2)
    
    # estimate with raoBust
    result_strong5 <- multinom_test(df$X, df$Y, strong=TRUE)
    
    # estimate with nnet 
    result_nnet5 <- nnet::multinom(cbind(df$Y[, 5], df$Y[, 1:4]) ~ 1)
    result_nnet5_summ <- summary(result_nnet5)
    
    # negative log likelihood 
    neg_ll_df[sim, 1] <- multinom_log_lik_alternative(as.vector(result_strong5$mle0), df$Y, df$X)
    neg_ll_df[sim, 2] <- multinom_log_lik_alternative(as.vector(t(result_nnet5_summ$coefficients)), df$Y, df$X)
    
    expect_equal(result_strong5$mle0[1, ], as.vector(result_nnet5_summ$coefficients), 
                 tolerance = 0.05)
  }
  
  summary(neg_ll_df$raoBust - neg_ll_df$nnet)
  
  # expect that nnet won't have a lower negative ll than raoBust every time
  expect_false(sum(neg_ll_df$raoBust - neg_ll_df$nnet < 0) == 0)
  
})

test_that("estimation under the weak null is not strictly worse than using the true B and matches the true B pretty well", {
  
  nsim <- 10
  nn <- 500
  
  neg_ll_df <- data.frame(raoBust = rep(NA, nsim),
                          true = rep(NA, nsim))
  
  for (sim in 1:nsim) {
    set.seed(sim)
    
    df <- simulate_data_mult(nn = nn, strong=FALSE, sd_beta1s=2, sd_beta0s=1, jj_null=2)
    
    # estimate with raoBust
    result_strong5 <- multinom_test(df$X, df$Y, strong=FALSE, j = 2)
    
    # negative log likelihood 
    neg_ll_df[sim, 1] <- multinom_log_lik_alternative(as.vector(result_strong5$mle0), df$Y, df$X)
    neg_ll_df[sim, 2] <- multinom_log_lik_alternative(as.vector(df$B), df$Y, df$X)
    
    expect_equal(result_strong5$mle0[1, ], df$B[1, ], 
                 tolerance = 0.05)
  }
  
  summary(neg_ll_df$raoBust - neg_ll_df$true)
  
  # expect that raoBust won't be much worse than the truth 
  expect_true(all.equal(neg_ll_df$raoBust, neg_ll_df$true, tolerance = 0.01*max(neg_ll_df$raoBust)))
  
})

