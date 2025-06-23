#' Test hypothesis of form A x beta = c
#'
#' @param test_object Object of type `raoFit`.
#'
#' @param A Matrix specifying linear combination of parameters for user-specified hypotheses.
#' 
#' @param c Vector that has user-specified null hypothesis value for user-specified linear combinations of parameters of interest.
#'
#' @return Table with relevant quantities of interest for hypothesis tests.
#'
#' @author Shirley Mathur
#'
#' @export

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