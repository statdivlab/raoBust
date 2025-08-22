set.seed(3)

### specifics not important, just for generating data
n_clust = 10; m = 3; sigma_b = 0.2; rho = 0.2
nn <- n_clust * m
F   <- rnorm(nn)
E   <- matrix(rnorm(nn * 4), nn, 4)
X   <- sqrt(rho) * F + sqrt(1 - rho) * E
X   <- scale(X, center = TRUE, scale = TRUE)
covariate1 <- X[, 1]; covariate2 <- X[, 2]
covariate3 <- X[, 3]; covariate4 <- X[, 4]
id  <- rep(1:n_clust, each = m)
b <- rnorm(n_clust, mean = 0, sd = sigma_b)

eta <- 1 + -1*covariate1 + 0*covariate2 + -1*covariate3 + 1*covariate4 + b[id]
mu  <- exp(eta)

yy  <- rpois(nn, lambda = mu)

df_re <- data.frame(
  id = id,
  covariate1, covariate2, covariate3, covariate4,
  yy = yy
)


expect_silent(glm_coef_default <- glm_test(
  formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
  family  = poisson(link = "log"),
  data    = df_re
)$coef)
expect_type(glm_coef_default, "double")


expect_silent(gee_coef_default <- gee_test(
  formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
  family  = poisson(link = "log"),
  id      = id,
  data    = df_re
)$coef)

expect_type(gee_coef_default, "list")

test_that("invariant to row order", {

  ### glm_test fine
  expect_equal(
    glm_coef_default,
    glm_test(
      formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      data    = df_re[sample(1:nrow(df_re), replace=FALSE), ]
    )$coef
  )

  expect_equal(
    gee_coef_default,
    gee_test(
      formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      id      = id,
      data    = df_re[nrow(df_re):1, ]
    )$coef
  )


  expect_equal(
    gee_coef_default,
    gee_test(
      formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      id      = id,
      data    = df_re[c(4:nrow(df_re), 1:3), ]
    )$coef
  )

  expect_equal(
    gee_coef_default,
    gee_test(
      formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      id      = id,
      data    = df_re[c(4:nrow(df_re), 3:1), ]
    )$coef
  )

  my_reorder <- sample(1:nrow(df_re), replace=FALSE)
  expect_equal(
    gee_coef_default,
    gee_test(
      formula = yy ~ covariate1 + covariate2 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      id      = id,
      data    = df_re[my_reorder, ]
    )$coef
  )

})


test_that("invariant to column order", {

  expect_equal(
    glm_coef_default,
    glm_test(
      formula = yy ~ covariate2 + covariate1 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      data    = df_re
    )$coef[c(1,3,2,4,5), ]
  )

  ### robust score test differs
  ###   FAILING
  expect_equal(
    gee_coef_default,
    gee_test(
      formula = yy ~ covariate2 + covariate1 + covariate3 + covariate4,
      family  = poisson(link = "log"),
      id      = id,
      data    = df_re
    )$coef[c(1,3,2,4,5,6), ]
  )
})
