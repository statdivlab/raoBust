#' Run robust wald test for null hypothesis of form A x beta = c
#'
#' @param test_object Object of type `raoFit`.
#'
#' @param A Matrix specifying linear combination of parameters for user-specified hypotheses.
#' 
#' @param c Vector that has user-specified null hypothesis value for user-specified linear combinations of parameters of interest.
#'
#' @return Table with relevant quantities of interest for hypothesis tests.
#' 
#' @importFrom stats rnorm
#'
#' @author Shirley Mathur
#'
#' @export
#' 
#' @examples
#' #set true value of beta for DGP
#' beta0s <- rnorm(n = 4)
#' beta1s <- rnorm(n = 4)
#' beta_true <- rbind(beta0s, beta1s)
#' beta_true[,2] <- 0 
#' beta_true[2,1] <- 0 
#' 
#' #generate sample data
#' sample_data <- simulate_data_mult(30, Beta = beta_true)
#' 
#' #run weak multinom test
#' sample_multinom_test <- multinom_test(X = sample_data$X,
#'                                       Y = sample_data$Y,
#'                                       j = 2)
#' 
#' #test hypothesis that the coefficient for covariate is same for j = 1, j = 2
#' #so, we test hypothesis that beta_{k=1,j=1} - beta_{k=1, j=2} = 0
#' 
#' #first, set up A matrix
#' my_A <- set_up_lin_com(J = 5, p = 1, n_hypotheses = 1) 
#' my_A[1,"k_1_j_1"] <- 1
#' my_A[1, "k_1_j_2"] <- -1
#' 
#' #then, set c
#' my_c <- 0
#'
#' #run linear combination test
#' sample_lincom_test <- lincom(sample_multinom_test,
#'                              A = my_A,
#'                              c = my_c)
#' 
#' #print result
#' sample_lincom_test
#' 

lincom <- function(test_object, A, c) {
  
  #get and compute relevant quantities
  n_hyp <- nrow(A)
  n <- nrow(test_object$design_matrix)
  c <- as.matrix(c, ncol = 1)
  beta_hat <- test_object$mle1
  wald_cov <- test_object$wald_cov
  beta_hat_vec <- matrix(as.vector(beta_hat), ncol = 1)
  s_beta <- A %*% beta_hat_vec - c
  
  #compute estimates under alternative
  lincom_est <- s_beta + c
  
  #compute std errors of estimates
  lincom_est_cov <- A %*% wald_cov %*% t(A)
  lincom_std_err <- sqrt(diag(lincom_est_cov))
  
  #compute upper and lower bounds for 95% confidence intervals for linear combinations
  lincom_upper <- lincom_est + lincom_std_err*qnorm(0.975)
  lincom_lower <- lincom_est - lincom_std_err*qnorm(0.975)
  
  #compute test statistic
  test_stat <- t(s_beta) %*% solve(lincom_est_cov) %*% s_beta

  #compute p-value
  pval <- 1 - pchisq(test_stat, nrow(A))
  
  #make coef table
  coef_tab <- data.frame("Est" = lincom_est,
                         "Std. Err" = lincom_std_err,
                         "Lower CI (95%)" = lincom_lower,
                         "Upper CI (95%)" = lincom_upper,
                         "pval" = pval)
  
  return(coef_tab)
  
}