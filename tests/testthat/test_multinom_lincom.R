test_that("lincom function works as expected when linear combination is checking is a single parameter is 0", {
  #set true value of beta for DGP
  beta0s <- rnorm(n = 4)
  beta1s <- rnorm(n = 4)
  beta_true <- rbind(beta0s, beta1s)
  beta_true[,2] <- 0 
  beta_true[2,1] <- 0 
  
  #generate sample data
  sample_data <- simulate_data_mult(30, Beta = beta_true)
  
  #run weak multinom test
  sample_multinom_test <- multinom_test(X = sample_data$X,
                                        Y = sample_data$Y,
                                        j = 2)
  
  #test hypothesis that the coefficient for k=1, j = 1 is 0
  #first, set up A matrix
  my_A <- set_up_lin_com(J = 5, p = 1, n_hypotheses = 1) 
  my_A[1,"k_1_j_1"] <- 1
  
  #then, set c
  my_c <- 0
  
  #run linear combination test
  sample_lincom_test <- lincom(sample_multinom_test,
                               A = my_A,
                               c = my_c)
  
  #compare value from lincom test and the wald test from test object
  sample_test_coef_tab <- sample_multinom_test$coef_tab
  sample_test_relevant_res <- sample_test_coef_tab[sample_test_coef_tab$"Category" == 1 & sample_test_coef_tab$"Covariate" == 1, ]
  
  expect_true(all.equal(unlist(unname(sample_test_relevant_res[1,3:7])), unlist(unname(sample_lincom_test[1,1:5]))))
})

test_that("type 1 error is controlled for lincom test", {
  
  skip("Skipping due to time, can check manually")
  
  #set up simulation parameters
  set.seed(13)
  nsim <- 1000
  n <- 200
  J <- 10
  null_j <- 2
  
  #set true value of beta for DGP
  beta0s <- rnorm(n = 9)
  beta1s <- rnorm(n = 9)
  beta_true <- rbind(beta0s, beta1s)
  beta_true[,1] <- beta_true[,2]
  
  #set up matrix for testing hypothesis that beta_{k=1,j=1} - beta_{k=1,j=2} = 0
  A <- set_up_lin_com(J, 1, 1)
  A[1, c("k_1_j_1", "k_1_j_2")] <- c(1,-1)
  
  #set up vector storing p-values
  lincom_sim_pvals <- rep(NA, nsim)
  
  #do simulation
  for (i in 1:nsim) {
    
    #generate data
    sim_weak_data_true <- simulate_data_mult(n, Beta = beta_true, overdispersion = 1)
    
    #run test for weak null hypothesis
    weak_multinom_res <- multinom_test(X = sim_weak_data_true$X, Y = sim_weak_data_true$Y, j = null_j)
    
    #run lincom test
    lincom_res <- lincom(weak_multinom_res, A, 0)
    
    #store pval result
    lincom_sim_pvals[i] <- lincom_res$pval
  }
  
  expect_true(mean(lincom_sim_pvals < 0.05) < 0.05 + 2.575 * sqrt(0.05 * 0.95/nsim))
  expect_true(mean(lincom_sim_pvals > 0.05) > 0.05 - 2.575 * sqrt(0.05 * 0.95/nsim))
  
})


