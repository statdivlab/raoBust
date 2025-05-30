set.seed(1)
data(cars)
cars2 <- cars
cars2$wts <- log(rnorm(nrow(cars), 20))
cars2$batch1 <- rep(LETTERS[1:5], each = 10)
cars2$batch2 <- rep(LETTERS[6:10], 10)
cars2$batch3 <- rep(1:5, each = 10)

test_that("replicates work", {

  ### fails in test unless `cars` is loaded?
  gee_test(formula = dist ~ speed,
           data = cars2,
           id = batch3,
           family=poisson(link="log"),
           offset = wts)

  ### Results should be consistent regardless of exact form of `id`
  ### Should be robust to ordering

  expect_identical(
    gee_test(formula = dist ~ speed,
             data = cars2,
             id = batch3,
             family=poisson(link="log"),
             offset = wts)$coef_tab,
    gee_test(formula = dist ~ speed,
             data = cars2,
             id = batch1,
             family=poisson(link="log"),
             offset = wts)$coef_tab
  )

  # devtools::load_all()
  expect_identical(
    gee_test(formula = dist ~ speed,
             data = cars2,
             id = batch2,
             family=poisson(link="log"),
             offset = wts)$coef_tab
    ,
    gee_test(formula = dist ~ speed,
             data = cars2[cars2$batch2 %>% order, ],
             id = batch2,
             family=poisson(link="log"),
             offset = wts)$coef_tab
  )

  expect_true({
    is.data.frame(gee_test(formula = dist ~ speed,
                           data = cars2,
                           id = batch1,
                           family=poisson(link="log"),
                           offset = wts)$coef_tab)
  })

  #### check gives warning if weights provided
  expect_warning(
    gee_test(formula = dist ~ speed,
             data = cars2,
             id = batch3,
             family=poisson(link="log"),
             offset = wts,
             weights = wts)
  )

  #### check correlation is estimated to be positive
  set.seed(2)
  cars3 <- rbind(cars, cars, cars)
  cars3$dist <- cars3$dist + sample(-2:2, replace=TRUE, size = length(cars3$dist))
  cars3$batch <- 1:length(cars$dist)
  expect_gt(gee_test(formula = dist ~ speed,
                       data = cars3,
                       id = batch,
                       family=poisson(link="log"))$coef_tab[3, "Estimate"], 0.9)


  #### check I've implemented score tests with clusters
  expect_false(
    any(is.na(gee_test(formula = dist ~ speed,
                       data = cars2,
                       id = batch3,
                       family=poisson(link="log"),
                       offset = wts)[["Robust Score p"]][1:2]))
  )
})

test_that("geeasy and geepack give similar results", {
  
  geeasy_test <- gee_test(formula = dist ~ speed,
                          data = cars2,
                          id = batch3,
                          family=poisson(link="log"),
                          offset = wts)$coef_tab
  
  geepack_test <- gee_test(use_geeasy = FALSE,
                           formula = dist ~ speed,
                           data = cars2,
                           id = batch3,
                           family=poisson(link="log"),
                           offset = wts)$coef_tab
  expect_true((geeasy_test$Estimate[1] - geepack_test$Estimate[1])/geeasy_test$Estimate[1] < 0.01)
})

test_that("jackknife standard errors work", {
  cars2_obj <- geelm(formula = dist ~ speed, data = cars2, id = batch3, 
                      family = poisson(link = "log"), offset = wts, corstr = "exchangeable")

  jack_se <- jackknife_se(cars2_obj, cars2)
  jack_se_cluster <- jackknife_se(cars2_obj, cars2, id = cars2$batch3)
  sand_se <- summary(cars2_obj)$coef$Std.err
  # expect that non-clustered se's are closer to robust se's because this data has no inherent clustering
  expect_true(sum((jack_se - sand_se)^2) < sum((jack_se_cluster - sand_se)^2))
  
  gee_res <- gee_test(formula = dist ~ speed,
                      data = cars2,
                      id = batch3,
                      family=poisson(link="log"),
                      offset = wts,
                      use_jack_se = TRUE)$coef_tab
  # expect that gee test robust standard errors are equal to those directly from jackknife function
  expect_true(all.equal(gee_res$`Robust Std Error`[1:2], as.vector(jack_se_cluster)))
})

test_that("glm solver works when gee solver fails", {
  gee_res <- gee_test(formula = dist ~ speed,
                      data = cars2,
                      id = batch3,
                      family=poisson(link="log"),
                      offset = wts)$coef_tab
  glm_res <- suppressWarnings(gee_test(formula = dist ~ speed,
                                       data = cars2,
                                       id = batch3,
                                       family=poisson(link="log"),
                                       offset = wts,
                                       skip_gee = TRUE))$coef_tab
  # check that estimates for coefficient for gee and glm result are similar 
  expect_equal(gee_res[2, 1], glm_res[2, 1], tolerance = 0.1)
  # confirm that without user-input cluster correlation coefficient, there
  # are no robust score p-values
  expect_true(is.na(glm_res[2, "Robust Score p"]))
  
  # check that warning is thrown when using glm instead of gee 
  expect_warning(glm_res_user_corr <- gee_test(formula = dist ~ speed,
                                               data = cars2,
                                               id = batch3,
                                               family=poisson(link="log"),
                                               offset = wts,
                                               skip_gee = TRUE,
                                               cluster_corr_coef = -0.04)$coef_tab)
  
  # check that robust score p-values are similar between gee and glm solver when
  # using gee estimated alpha parameter as user-input cluster correlation coefficient
  expect_equal(gee_res[2, "Robust Score p"], glm_res_user_corr[2, "Robust Score p"], tolerance = 0.01)
})
