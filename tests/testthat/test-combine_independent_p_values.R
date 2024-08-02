test_that("uniform p-values under H0", {
  xx <- replicate(50000,
                  {
                    combine_independent_p_values(runif(3))
                  })
  expect_equal(mean(xx < 0.1), 0.1, tolerance=0.02)
  expect_equal(mean(xx < 0.5), 0.5, tolerance=0.02)
  expect_equal(mean(xx < 0.9), 0.9, tolerance=0.02)
# hist(xx)
})

test_that("non-trivial power", {
  ## simulate small p-values only, between 0 and 0.1
  xx <- replicate(50000,
                  {
                    combine_independent_p_values(runif(3, 0, 0.10))
                  })
  expect_true(mean(xx < 0.1) > 0.8)
  expect_true(mean(xx < 0.01) > 0.5)
  # hist(xx)
})
