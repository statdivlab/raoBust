#' Negative log-likelihood for multinomial data under the weak null (\eqn{\beta_j = 0} for specific j)
#'
#'
#' @param beta_as_vector A vector containing the values for all \eqn{\beta_k}, with \eqn{k \neq j}, as well as all \eqn{\beta_{k0}, for k = 1, \dots, J}.
#' In particular, this vector should be so that the first \eqn{(J-2)(p+1)} entries are \eqn{\beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots, \beta_{(j-1)0}, \beta_{j-1}^{\top},  \beta_{(j+1)0}, \beta_{j+1}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}}.
#' Then, the \eqn{(J-2)(p+1) + 1} entry should be \eqn{\beta_{j0}}.
#' @param beta_j_null This is the null hypothesized value for \eqn{\beta_j}, which is by default set to be \eqn{\beta_j = 0}.
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param j This specifies for which category you want to compute the MLE under the constraint that \eqn{\beta_j = 0}.
#'
#' @return The value of the log likelihood for the input \eqn{\beta}
#'
#' @author Shirley Mathur
#'
#' @export
multinom_log_lik_weak_null <- function(beta_as_vector, Y, X, j, beta_j_null = NULL) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  p <- ncol(X) #get number of covariates
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  if (is.null(beta_j_null)) {
    beta_j_null <- rep(0, p)
  }

  # beta <- beta_vector_to_matrix(values = beta_as_vector, p = p, J = J) ## (p+1) x (J-1) beta matrix
  # beta[2:(p+1), j] <- beta_j_null
  
  beta <- beta_vector_to_matrix(beta_as_vector, p, J, j, beta_j_null) 
  
  # beta <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
  # for (k in 1:(J-1)) {
  #   if (k < j) {
  #     beta[ ,k] <- beta_as_vector[((k-1)*(p+1)+1):(k*(p+1))]
  #   } else if (k == j) {
  #     beta[,k] <- c(beta_as_vector[(p+1)*(J-2)+1], beta_j_null)
  #   } else if (k > j) {
  #     beta[ ,k] <- beta_as_vector[((k-2)*(p+1)+1):((k-1)*(p+1))]
  #   }
  # }

  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  XaugBeta <- Xaug %*% beta
  pJ <- (1 + rowSums(exp(XaugBeta)))^(-1)
  ps <- pJ*exp(XaugBeta)
  ps_full <- cbind(ps, pJ)

  if (any(ps_full <= 0)) {
    # print("beta")
    # print(beta)
    # print("ps_full")
    # print(ps_full)
    # stop("Some probabilities are estimated to be zero")
  }


  # if (any(ps_full > 0)) {
  #   loglik <- .Machine$double.xmax
  # } else {
  loglik <- sum(Y*log(ps_full))
  # }

  return(-loglik)
}


