nn <- 100

simulate_null_data <- function() {
  covariate1 <- rnorm(nn)
  batch_dev <- data.frame("batch" = 1:(nn/5), "deviation" = rnorm(nn/5, 0, 1))
  means <- merge(expand.grid("obs_index" = 1:5, "batch" = 1:(nn/5)), batch_dev, by = "batch")
  means$mean <- 1 + 0*covariate1 + means$deviation
  yy <- rpois(nn, lambda = exp(means$mean))
  df <- data.frame(covariate1 = covariate1,
                   batch = means$batch,
                   yy = yy)
  return(df)
}
test_that("nominal level generalized score test for n = 100 under correct model specification", {

  set.seed(2)
  nsim <- 200
  ps <- vector("numeric", length = nsim)

  expect_silent({
    for (j in 1:nsim) {
      ps[j] <- gee_test(formula = yy ~ covariate1,
                        data = simulate_null_data(),
                        id = batch,
                        family = poisson(link = "log"))$coef_tab[2, "Robust Score p"]
    }
  })

  # for (j in 1:nsim) {
  #   ps[j] <- glm_test(yy ~ covariate1,
  #                     data = simulate_null_data(),
  #                     family = poisson(link = "log"))[2, "Robust Score p"]
  # }

  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10))

  ### for 100 batches of 5, (1%, 5%, 10%) error is (0.8, 4.8, 10.8)% with no correction (500 sims)
  ### for 20 batches of 5, it's                    (0.0, 5.2, 13.0)% with no correction (500 sims)

  expect_true(mean(ps < 0.05) < 0.05 + 2.575 * sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.05) > 0.05 - 2.575 * sqrt(0.05 * 0.95/nsim))

})
