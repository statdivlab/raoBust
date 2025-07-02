#' Optimization under null or alternative for multinomial model via Fisher scoring.
#'
#' @param beta The initial values provided for the \eqn{\beta} parameters.
#' @param Y The \eqn{n x J} data matrix of outcomes.
#' @param X The \eqn{n x p} design matrix of covariates.
#' @param null If TRUE, optimizes under the null, if FALSE, optimizes under the alternative. Defaults to TRUE.
#' @param strong If FALSE, this function will compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0} for the length \eqn{p} vector \eqn{\beta_j}.
#' If TRUE, this function instead computes the robust score statistic to test the strong null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0} for all length \eqn{p} vectors \eqn{\beta_j}, \eqn{j\in\{1,\ldots,J-1\}}. 
#' Default is FALSE.
#' @param null_j If `strong` is FALSE, this argument must be supplied. This gives the category \eqn{j} in the weak null hypothesis that \eqn{\beta_j = 0}.
#' @param tol The tolerance used to determine how much better update function value must be prior to stopping algorithm.
#' @param stepSize The size of the step to take during the parameter update step.
#' @param arm_c Control parameter for checking Armijo condition.
#' @param maxit Maximum number of iterations for Fisher scoring. Defaults to 250.
#' @param pseudo_inv Use the pseudo-inverse of the Fisher information matrix for the update (in case the inverse in computationally singular)
#' @return The optimal beta values under the null or alternative model.
#'
#' @author Shirley Mathur
#'
#' @export
multinom_fisher_scoring <- function(beta, X, Y, null = TRUE, strong = FALSE, null_j = NULL, tol = 1e-5, stepSize = 0.5, arm_c = 0.5, maxit = 250, pseudo_inv = FALSE) {
  
  # check that if strong is FALSE, j is provided 
  if (null & !strong) {
    if (is.null(null_j)) {
      stop("If testing under the weak null hypothesis with `strong = FALSE`, you must include the argument `null_j`.")
    }
  }
  
  #get relevant quantities needed for computations
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)
  
  #compute other necessary quantities needed for computations
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  
  
  #get indices of beta components that need to be optimized over
  if (null) {
    if (strong == FALSE) {
      optim_indices <-  (1:((p+1)*(J-1)))[-c(((null_j-1)*(p+1)+2):(null_j*(p+1)))] #take all indices that are not set to 0 under weak null
      optim_indices <- c(optim_indices[-((null_j-1)*(p+1) + 1)], optim_indices[((null_j-1)*(p+1) + 1)])
    } else {
      optim_indices <- (1:((p+1)*(J-1)))[((0:(J-2))*(p+1)+1)] #take all indices that are not set to 0 under strong null
    }
  } else {
    optim_indices <- (1:((p+1)*(J-1)))
  }
  
  #initialize beta values, beta matrix, probs matrix, and likelihood for initial beta value
  beta_current <- beta
  beta_values_current <- as.vector(beta_current)[optim_indices]
  probs_current <- multinom_get_probs(X, Y, beta)
  logLik_current <- sum(Y*log(probs_current))

  #initialize converged as FALSE
  converged <- FALSE
  n_it <- 1

  #do algorithm for finding optimal values of beta
  while (!converged & n_it < maxit) {
    #compute score vector
    score <- multinom_score_vector(X, Y, probs_current)
    score <- score[optim_indices]
    
    #compute information matrix
    info_mat <- multinom_info_mat(X, Y, probs_current)
    info_mat <- info_mat[optim_indices, optim_indices]
    
    #compute step direction
    if (!pseudo_inv) {
      step_dir <- tryCatch({solve(info_mat) %*% score},
                           error = function(cond) {
                             print(cont)
                             return(NA)
                            })
    } else {
      step_dir <- tryCatch({solve(info_mat) %*% score},
                           error = function(cond) {
                             return(MASS::ginv(info_mat) %*% score)
                           })
    }
    
    #first get acceptable update
    accepted <- FALSE
    
    while(!accepted) {
      #compute proposed beta
      beta_values_prop <- beta_values_current + stepSize*step_dir

      if (null) {
        if (strong == FALSE) {
          beta_prop <- multinom_beta_vector_to_matrix(beta_values_prop, p, J, null_j)
        } else {
          beta_prop <- matrix(0, nrow = p + 1, ncol = J - 1)
          beta_prop[1,] <- beta_values_prop
        }
      } else {
        beta_prop <- matrix(beta_values_prop, nrow = p + 1, ncol = J - 1)
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
    
    n_it <- n_it + 1
    if (n_it > maxit) {
      message("Maximum number of iterations for MLE reached, exiting optimization")
    }
  }
  
  beta_null_mle <- beta_current
  
  return(beta_null_mle)
    
}
