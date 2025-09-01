test_that("test Poisson example", {

  expect_type(glm_test(dist ~ speed, data = cars, family=poisson(link="log"))$coef_tab,
              "double")

  expect_error(glm_test(dist ~ speed, data = cars))


  expect_type(glm_test(dist ~ speed,
                       data = cbind(cars, "offs" = rnorm(nrow(cars), 20)),
                       family=poisson(link="log"),
                       offset = log(offs))$coef_tab,
              "double")

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"))$coef_tab,
              "double")

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20)))$coef_tab,
              "double")
  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20)),
                       weights = NULL)$coef_tab,
              "double")

  expect_error(glm_test(dist ~ speed,
                        data = cars,
                        family=poisson(link="log"),
                        offset = log(rnorm(nrow(cars), 20)),
                        weights = rep(1, 20)))

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20)))$coef_tab,
              "double")

})

test_that("test logistic example", {
  
  expect_type(glm_test((dist > 43) ~ speed, data = cars, family=binomial(link = "logit"))$coef_tab,
              "double")
  
  cars$id <- rep(1:5, 10)
  expect_type(gee_test(formula = (dist > 43) ~ speed, data = cars, id = "id", 
                       family=gaussian(link = "identity"))$coef_tab,
              "list")
  
})

test_that("test gaussian example", {
  
  expect_type(glm_test(dist ~ speed, data = cars, family=gaussian(link = "identity"))$coef_tab,
              "double")
  
  cars$id <- rep(1:5, 10)
  expect_type(gee_test(formula = dist ~ speed, data = cars, id = "id", 
                       family=gaussian(link = "identity"))$coef_tab,
              "list")
  
})

test_that("error with other model", {
  
  expect_error(glm_test(dist ~ speed, data = cars, family=gaussian(link = "log")))

})

