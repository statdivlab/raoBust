#' Optimization under null constraint for multinomial model via Fisher scoring.
#'
#'@param beta This should be the initial values provided for the \eqn{\beta} parameters.
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param strong This is by default specified as FALSE to compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0}.
#' If specified to be TRUE, the function instead computes the robust score statistic to test the strong null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0}.
#' @param j If `strong` is specified as FALSE, this argument must be supplied. This specifies for which category \eqn{j} you want to test the weak null hypothesis that \eqn{\beta_j = 0}.
#' @param tol This is tolerance used to determine how much better update function value must be prior to stopping algorithm.
#' @param stepSize This is the size of the step to take during the parameter update step.
#' @param c Control parameter for checking Armijo condition.
#' @return The optimal beta values under the null constraint.
#'
#'
#' @author Shirley Mathur
#'
#' @export
multinom_fisher_scoring <- function(beta, X, Y, strong = FALSE, null_j = NULL, tol = 1e-5, stepSize = 0.5, arm_c = 0.5) {
  
  
  #get relevant quantities needed for computations
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)
  
  #compute other necessary quantities needed for computations
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  
  
  #get indices of beta components that need to be optimized over
  if (strong == FALSE) {
    optim_indices <-  (1:((p+1)*(J-1)))[-c(((null_j-1)*(p+1)+2):(null_j*(p+1)))] #take all indices that are not set to 0 under weak null
    optim_indices <- c(optim_indices[-((null_j-1)*(p+1) + 1)], optim_indices[((null_j-1)*(p+1) + 1)])
  } else {
    optim_indices <- (1:((p+1)*(J-1)))[((0:(J-2))*(p+1)+1)] #take all indices that are not set to 0 under strong null
  }
  
  #initialize beta values, beta matrix, probs matrix, and likelihood for initial beta value
  beta_current <- beta
  beta_values_current <- as.vector(beta_current)[optim_indices]
  probs_current <- multinom_get_probs(X, Y, beta)
  logLik_current <- sum(Y*log(probs_current))

  #initialize converged as FALSE
  converged <- FALSE

  #do algorithm for finding optimal values of beta
  while (!converged) {
    #compute score vector
    score <- multinom_score_vector(X, Y, probs_current)
    score <- score[optim_indices]
    
    #compute information matrix
    info_mat <- multinom_info_mat(X, Y, probs_current)
    info_mat <- info_mat[optim_indices, optim_indices]
    
    #compute step direction
    step_dir <- solve(info_mat) %*% score
    
    #first get acceptable update
    accepted <- FALSE
    
    while(!accepted) {
      #compute proposed beta
      beta_values_prop <- beta_values_current + stepSize*step_dir

      if (strong == FALSE) {
        beta_prop <- beta_vector_to_matrix(beta_values_prop, p, J, null_j)
      } else {
        beta_prop <- matrix(0, nrow = p + 1, ncol = J - 1)
        beta_prop[1,] <- beta_values_prop
      }
      
      #compute likelihood for proposed beta
      probs_prop <- multinom_get_probs(X, Y, beta_prop)
      logLik_prop <- sum(Y*log(probs_prop))
      
      #check if update is acceptable
      if (logLik_prop >= logLik_current + stepSize*arm_c*norm(score %*% step_dir)) {
        beta_values_update <- beta_values_prop
        beta_update <- beta_prop
        probs_update <- probs_prop
        logLik_update <- logLik_prop
        accepted <- TRUE
      } else {
        stepSize <- 0.5*stepSize
      }
    }
    
    #once we get acceptable update, check if converged
    rel_diff <- abs((logLik_current - logLik_update)/(logLik_current))
    
    if (rel_diff < tol) {
      beta_null_mle <- beta_current
      converged <- TRUE
    } else {
      beta_current <- beta_update
      beta_values_current <- beta_values_update
      probs_current <- probs_update
      logLik_current <- logLik_update
    }
  }
  
  beta_null_mle <- beta_current
  
  return(beta_null_mle)
    
}
