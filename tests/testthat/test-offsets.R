my_dff <- data.frame("x0" = 1:50,
                    "offs" = 1/(1:50))
my_dff$y0 <- pmax(0.1, rnorm(50, mean = 10*my_dff$x0 * my_dff$offs, sd = 1))

test_that("offsets impact the test and estimates", {

  #### expect close to 1 and zero

  expect_silent({
    glm_off <- glm_test(y0 ~ log(x0),
                        data = my_dff,
                        family=poisson(link="log"),
                        offset = log(offs))$coef_tab
    glm_no_off <- glm_test(y0 ~ log(x0),
                           data = my_dff,
                           family=poisson(link="log"))$coef_tab
  })

  expect_equal(glm_off[2, 1],
               expected= 1,
               tolerance=0.05)

  expect_equal(glm_no_off[2, 1],
               expected=0,
               tolerance=0.05)

  expect_true(glm_off[2, "Robust Score p"] != glm_no_off[2, "Robust Score p"])

})
