test_that("see what happens with missing data", {
  
  set.seed(1)
  n_clust = 10
  m = 3
  sigma_b = 0.1
  rho = 0.2
  nn <- n_clust * m
  
  #### generate covariates to be correlated. rho = 0 means uncorrelated.
  F   <- rnorm(nn)
  E   <- matrix(rnorm(nn * 4), nn, 4)
  X   <- sqrt(rho) * F + sqrt(1 - rho) * E
  X   <- scale(X, center = TRUE, scale = TRUE)
  
  covariate1 <- X[, 1]; covariate2 <- X[, 2]
  covariate3 <- X[, 3]; covariate4 <- X[, 4]
  
  #### generate observations to be cluster correlated via random effect
  id  <- rep(1:n_clust, each = m)
  b <- rnorm(n_clust, mean = 0, sd = sigma_b)
  
  #### can toggle effect sizes here
  eta <- 1 + -3*covariate1 - 0*covariate2 + -1*covariate3 + 1*covariate4 + b[id]
  mu  <- exp(eta)
  
  yy  <- rpois(nn, lambda = mu)
  
  yy[5] <- NA # TODO stress test
  
  df_re <- data.frame(
    id = id,
    covariate1, covariate2, covariate3, covariate4,
    yy = yy
  )
  
  # see what happens for glm when missing y value
  res_all <- glm_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    data    = df_re
  )
  res_no_5 <- glm_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    data    = df_re[-5,]
  )
  expect_true(all.equal(res_all$coef_tab, res_no_5$coef_tab))
  
  # how about gee? 
  gee_res_all <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re
  )
  gee_res_no5 <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re[-5, ]
  )
  gee_res_geelm <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re,
    use_geeasy = FALSE
  )
  gee_res_geelm_no5 <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re[-5, ],
    use_geeasy = FALSE
  )
  all.equal(gee_res_no5$coef_tab, gee_res_all$coef_tab)
  all.equal(gee_res_geelm$coef_tab, gee_res_geelm_no5$coef_tab)
  
  # what if it is covariate that is missing? 
  yy  <- rpois(nn, lambda = mu)
  
  covariate2[4] <- NA # TODO stress test
  
  df_re <- data.frame(
    id = id,
    covariate1, covariate2, covariate3, covariate4,
    yy = yy
  )
  
  # see what happens for glm when missing y value
  res_all <- glm_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    data    = df_re
  )
  res_no_4 <- glm_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    data    = df_re[-4,]
  )
  expect_true(all.equal(res_all$coef_tab, res_no_4$coef_tab))
  
  # how about gee? 
  gee_res_all <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re
  )
  gee_res_no4 <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re[-4, ]
  )
  all.equal(gee_res_no4$coef_tab, gee_res_all$coef_tab)
  gee_res_geelm <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re,
    use_geeasy = FALSE
  )
  gee_res_no4_geelm <- gee_test(
    formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
    family  = poisson(link = "log"),
    id      = id,
    data    = df_re[-4, ],
    use_geeasy = FALSE
  )
  all.equal(gee_res_no4_geelm$coef_tab, gee_res_geelm$coef_tab)
})
