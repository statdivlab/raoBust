
nsim <- 200
nn <- 100
simulate_null_data <- function() {
  covariate1 <- rnorm(nn)
  yy <- rpois(nn, lambda = 5 + 0*covariate1)
  df <- data.frame(covariate1 = covariate1,
                   yy = yy)
  return(df)
}

test_that("nominal level score test for n = 100 under correct model specification", {

  set.seed(1)
  nsim <- 200
  ps <- vector("numeric", length = nsim)

  expect_silent({
    for (j in 1:nsim) {
      ps[j] <- glm_test(yy ~ covariate1,
                        data = simulate_null_data(),
                        family = poisson(link = "log"))$coef_tab[2, "Robust Score p"]
    }
  })

  expect_true(mean(ps < 0.05) < 0.05 + 2.575 * sqrt(0.05 * 0.95/nsim))
  expect_true(mean(ps < 0.05) > 0.05 - 2.575 * sqrt(0.05 * 0.95/nsim))

})

test_that("p-values are correlated with robust Wald under null", {

  set.seed(2)
  nsim <- 200
  ps_score <- vector("numeric", length=nsim)
  ps_wald <- vector("numeric", length=nsim)
  nn <- 100

  for (i in 1:nsim) {
    the_glm <- glm_test(yy ~ covariate1,
                        data = simulate_null_data(),
                        family = poisson(link = "log"))$coef_tab
    ps_score[i] <- the_glm[2, "Robust Score p"]
    ps_wald[i] <- the_glm[2, "Robust Wald p"]
  }

  lm_res <- coef(summary(lm(ps_score ~ ps_wald)))

  expect_equal(lm_res[1,1], 0, tolerance = 0.02)
  expect_equal(lm_res[2,1], 1, tolerance = 0.02)

  expect_true(cor(ps_score, ps_wald) > 0.97)

})

simulate_alt_data <- function(nn) {
  covariate1 <- rnorm(nn)
  yy <- rpois(nn, lambda = 5 + 0.2*covariate1)
  df <- data.frame(covariate1 = covariate1,
                   yy = yy)
  return(df)
}

test_that("score test has increasing power with n = 100 under correct model specification", {

  nsim <- 200

  set.seed(1)
  nn <- 50
  ps <- vector("numeric", length = nsim)

  expect_silent({
    for (j in 1:nsim) {
      ps[j] <- glm_test(yy ~ covariate1,
                        data = simulate_alt_data(nn),
                        family = poisson(link = "log"))$coef_tab[2, "Robust Score p"]
    }
  })
  power50 <- mean(ps < 0.05)

  nn <- 100
  ps <- vector("numeric", length = nsim)

  expect_silent({
    for (j in 1:nsim) {
      ps[j] <- glm_test(yy ~ covariate1,
                        data = simulate_alt_data(nn),
                        family = poisson(link = "log"))$coef_tab[2, "Robust Score p"]
    }
  })
  power100 <- mean(ps < 0.05)

  nn <- 200
  ps <- vector("numeric", length = nsim)

  expect_silent({
    for (j in 1:nsim) {
      ps[j] <- glm_test(yy ~ covariate1,
                        data = simulate_alt_data(nn),
                        family = poisson(link = "log"))$coef_tab[2, "Robust Score p"]
    }
  })
  power200 <- mean(ps < 0.05)

  expect_true(power50 <= power100)
  expect_true(power100 <= power200)

})
