age <- c(72, 81, 90, 72, 70, 72, 75, 75, 67, 70, 86, 71, 75, 70, 79,
         70, 81, 72, 68, 67, 84, 71, 70, 71, 89, 78, 71, 86, 71, 70, 71,
         70, 73, 73, 67, 80, 74, 75, 70, 72, 81, 70, 71, 70, 69, 82, 74,
         74, 72, 68, 73, 84, 79, 76, 86, 71, 67, 75, 72, 75, 78, 70, 74,
         68, 73, 70, 74, 86, 70, 72, 72, 76, 68, 70, 73, 79, 72, 68, 71,
         67, 71, 68, 81, 83, 73, 79, 69, 82, 77, 84, 81, 83, 80, 76, 75,
         79, 75, 72, 75, 79)
sex <- c("Male", "Female", "Male", "Male", "Female", "Male", "Male",
         "Male", "Female", "Female", "Female", "Male", "Female", "Female",
         "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Male",
         "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male",
         "Male", "Female", "Male", "Female", "Male", "Female", "Female",
         "Female", "Female", "Female", "Female", "Female", "Female", "Male",
         "Male", "Male", "Female", "Male", "Male", "Male", "Female", "Male",
         "Male", "Male", "Female", "Male", "Male", "Female", "Female",
         "Male", "Female", "Female", "Male", "Male", "Male", "Female",
         "Female", "Male", "Male", "Male", "Female", "Female", "Female",
         "Male", "Female", "Male", "Female", "Male", "Male", "Female",
         "Male", "Female", "Male", "Female", "Male", "Male", "Male", "Male",
         "Female", "Male", "Female", "Male", "Female", "Female", "Male",
         "Female", "Male", "Female", "Male", "Female")
weight <- c(173, 139, 145, 190, 153, 154.5, 161.5, 158.2, 168, 127, 139,
            181.2, 137, 137, 181, 180.5, 151.5, 171, 210, 160, 154, 191,
            181, 124, 111, 172.5, 167, 167, 152.5, 183.5, 204, 135, 208,
            168.5, 230, 88, 163, 86, 167, 167, 158, 192, 126, 155, 239, 146,
            122, 162, 204.2, 168, 136, 156, 170, 193.5, 182, 155, 196, 145.5,
            135, 159, 118, 141, 241, 158.5, 152, 107, 171, 156, 177, 152,
            166, 159, 137, 131, 145, 174, 170, 135, 170.5, 142, 186, 217,
            161, 98, 124, 176, 198, 180.5, 152, 128, 137, 153.5, 130, 97,
            163, 179, 150, 157, 164.5, 139)

mri_mini <- data.frame(age, sex, weight)

test_that("returns a real test statistic", {

  r1 <- glm(age ~ sex + weight, data = mri_mini, family=poisson(link="log")) # it doesn't need to make sense

  pp1 <- robust_score_test(age ~ sex + weight, df = mri_mini, 1)
  pp2 <- robust_score_test(age ~ sex + weight, df = mri_mini, 2)

  expect_true(min(pp1, pp2) > 0)
  expect_true(max(pp1, pp2) < 1)

})

test_that("test multiple parameters at once", {

  pp3 <- robust_score_test(age ~ sex + weight, df = mri_mini, 1:2)

  expect_true(pp3 > 0)
  expect_true(pp3 < 1)

})

test_that("nominal level test for n=100 under correct model specification", {

  set.seed(1)
  nsim = 500
  ps <- vector("numeric", length=nsim)
  nn <- 100

  for (i in 1:nsim) {
    covariate1 <- rnorm(nn)
    yy <- rpois(nn, lambda = 5 + 0*covariate1)
    ps[i] <- robust_score_test(yy ~ covariate1, df = data.frame(covariate1, yy), 2)

  }
  expect_true(mean(ps < 0.05) < 0.05 + 2.575*sqrt(0.05*0.95/nsim))
  expect_true(mean(ps < 0.05) > 0.05 - 2.575*sqrt(0.05*0.95/nsim))


})

test_that("p-values are correlated with robust Wald under null", {

  set.seed(2)
  nsim = 200
  ps_score <- vector("numeric", length=nsim)
  # ps_rwald <- vector("numeric", length=nsim)
  ps_wald <- vector("numeric", length=nsim)
  nn <- 100

  for (i in 1:nsim) {
    covariate1 <- rnorm(nn)
    yy <- rpois(nn, lambda = 5 + 0*covariate1)
    ps_score[i] <- robust_score_test(yy ~ covariate1, df = data.frame(covariate1, yy), 2)
    # ps_rwald[i] <- coef(regress("rate", yy ~ covariate1, data = data.frame(covariate1, yy)))[2, "Pr(>|t|)"]
    ps_wald[i] <- coef(summary(glm(yy ~ covariate1, data = data.frame(covariate1, yy)), family=poisson(link="log")))[2, "Pr(>|t|)"]

  }

  lm_res <- coef(summary(lm(ps_score ~ ps_wald)))
  expect_equal(lm_res[1,1], 0, tolerance=0.02)
  expect_equal(lm_res[2,1], 1, tolerance=0.02)


  # expect_true(cor(ps_score, ps_rwald) > 0.97)
  expect_true(cor(ps_score, ps_wald) > 0.97)

})
