test_that("offsets work", {

  expect_type(glm_test(dist ~ speed, data = cars, family=poisson(link="log")),
              "double")

  expect_error(glm_test(dist ~ speed, data = cars))


  expect_type(glm_test(dist ~ speed,
                       data = cbind(cars, "offs" = rnorm(nrow(cars), 20)),
                       family=poisson(link="log"),
                       offset = log(offs)),
              "double")

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log")),
              "double")

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20))),
              "double")
  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20)),
                       weights = NULL),
              "double")

  expect_error(glm_test(dist ~ speed,
                        data = cars,
                        family=poisson(link="log"),
                        offset = log(rnorm(nrow(cars), 20)),
                        weights = rep(1, 20)))

  expect_type(glm_test(dist ~ speed,
                       data = cars,
                       family=poisson(link="log"),
                       offset = log(rnorm(nrow(cars), 20))),
              "double")

})
