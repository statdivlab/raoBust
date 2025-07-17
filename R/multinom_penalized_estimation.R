#' Optimization under null or alternative for multinomial model with a Firth penalty.
#'
#' @param beta The initial values provided for the \eqn{\beta} parameters.
#' @param Y The \eqn{n x J} data matrix of outcomes.
#' @param X The \eqn{n x p} design matrix of covariates.
#' @param null If TRUE, optimizes under the null, if FALSE, optimizes under the alternative. Defauls to TRUE.
#' @param strong If FALSE, this function will compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0} for the length \eqn{p} vector \eqn{\beta_j}.
#' If TRUE, this function instead computes the robust score statistic to test the strong null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0} for all length \eqn{p} vectors \eqn{\beta_j}, \eqn{j\in\{1,\ldots,J-1\}}. 
#' Default is FALSE.
#' @param null_j If `strong` is FALSE, this argument must be supplied. This gives the category \eqn{j} in the weak null hypothesis that \eqn{\beta_j = 0}.
#' @param j_ind If `strong` is FALSE and `null_j` is NULL, this argument must be supplied. This gives the category index of the individual covariate that is tested in the weak null hypothesis that \eqn{\beta_{kj} = 0}.
#' @param k_ind If `strong` is FALSE and `null_j` is NULL, this argument must be supplied. This gives the covariate index of the individual covariate that is tested in the weak null hypothesis that \eqn{\beta_{kj} = 0}.
#' @param tol The tolerance used to determine how much better update function value must be prior to stopping algorithm.
#' @param stepSize The size of the step to take during the parameter update step, used in the MLE step.
#' @param arm_c Control parameter for checking Armijo condition, used in the MLE step.
#' @param maxit Maximum number of iterations for augmentation algorithm. Defaults to 250.
#' @param maxit_fs Maximum number of iterations for Fisher scoring. Defaults to 5.
#' @param pseudo_inv Use the pseudo-inverse of the Fisher information matrix for the update (in case the inverse in computationally singular)
#' @return The optimal beta values under the null or alternative model.
#'
#' @author Sarah Teichman
#'
#' @export
multinom_penalized_estimation <- function(beta, X, Y, null = TRUE, strong = FALSE, null_j = NULL, j_ind = NULL, k_ind = NULL, 
                                          tol = 1e-5, stepSize = 0.5, arm_c = 0.5,
                                          maxit = 250, maxit_fs = 5, pseudo_inv = FALSE) {
  
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
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  
  # calculate H matrix for augmentations
  # get expanded X matrix
  X_cup <- multinom_X_cup_from_X(Xaug, J)
  
  # get G matrix
  G <- multinom_get_G(Xaug, J, X_cup)
  
  #initialize converged as FALSE and n_it as 1
  converged <- FALSE
  n_it <- 1
  
  Y_curr <- Y
  beta_curr <- beta
  
  while (!converged & n_it < maxit) {
    
    # set new z values
    z <- log(rowSums(Y)) - log(rowSums(exp(Xaug %*% beta_curr)))
    
    # calculate augmentations
    augs <- multinom_get_augmentations(G, Y, beta_curr, z, pseudo_inv)
    Y_curr <- Y + augs
    
    # get new beta values
    beta_old <- beta_curr
    beta_curr <- multinom_fisher_scoring(beta = beta_old, X = X, Y = Y_curr,
                                         null = null, strong = strong, null_j = null_j, j_ind = j_ind, k_ind = k_ind,
                                         tol = tol, stepSize = stepSize, arm_c = arm_c,
                                         maxit = maxit_fs, pseudo_inv = pseudo_inv)
    
    # check for convergence 
    B_diff <- max(abs(beta_curr - beta_old))
    if (B_diff < tol) {
      converged <- TRUE
    }
    
    n_it <- n_it + 1
    
    if (n_it > maxit) {
      message("Maximum number of iterations has been reached for penalty augmentation algorithm, exiting optimization.")
    }
    
  }
  
  return(beta_curr)
  
}
