test_that("replicates work", {

  set.seed(1)
  cars2 <- cars
  cars2$wts <- log(rnorm(nrow(cars), 20))
  cars2$batch1 <- rep(LETTERS[1:5], each = 10)
  cars2$batch2 <- rep(LETTERS[6:10], 10)
  cars2$batch3 <- rep(1:5, each = 10)

  gee_test(dist ~ speed,
           data = cars2,
           id = batch3,
           family=poisson(link="log"),
           offset = wts)

  ### Results should be consistent regardless of exact form of `id`
  ### Should be robust to ordering

  expect_identical(
    gee_test(dist ~ speed,
             data = cars2,
             id = batch3,
             family=poisson(link="log"),
             offset = wts),
    gee_test(dist ~ speed,
             data = cars2,
             id = batch1,
             family=poisson(link="log"),
             offset = wts)
  )

  # devtools::load_all()
  expect_identical(
    gee_test(dist ~ speed,
             data = cars2,
             id = batch2,
             family=poisson(link="log"),
             offset = wts)
    ,
    gee_test(dist ~ speed,
             data = cars2[cars2$batch2 %>% order, ],
             id = batch2,
             family=poisson(link="log"),
             offset = wts)
  )

  expect_true({
    is.data.frame(gee_test(dist ~ speed,
                           data = cars2,
                           id = batch1,
                           family=poisson(link="log"),
                           offset = wts))
  })

  #### check gives warning if weights provided
  expect_warning(
    gee_test(dist ~ speed,
             data = cars2,
             id = batch3,
             family=poisson(link="log"),
             offset = wts,
             weights = wts)
  )

  #### check correlation is estimated to be positive
  expect_true(FALSE)

  #### I've implemented score tests
  expect_false(
    any(is.na(gee_test(dist ~ speed,
                       data = cars2,
                       id = batch3,
                       family=poisson(link="log"),
                       offset = wts)[["Robust Score p"]]))
  )
})

