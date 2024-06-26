#' Negative log-likelihood for multinomial data under the alternative
#'
#'
#' @param beta_as_vector A vector containing the values for all \eqn{\beta_k} for \eqn{k = 1, \dots, J-1} and \eqn{\beta_{k0}, for k = 1, \dots, J}.
#' In particular, this vector should be so that the \eqn{(J-1)(p+1)} entries are \eqn{\beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}}.
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#'
#' @return The value of the log likelihood for the input \eqn{\beta}
#'
#' @author Shirley Mathur
#'
#' @export
multinom_log_lik_alternative <- function(beta_as_vector, Y, X) {
  
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  p <- ncol(X) #get number of covariates
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  
  
  beta <- matrix(beta_as_vector, nrow = p+1, ncol = J-1)  #initialize full (p+1) x (J-1) beta matrix
  
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  XaugBeta <- Xaug %*% beta
  pJ <- (1 + rowSums(exp(XaugBeta)))^(-1)
  ps <- pJ*exp(XaugBeta)
  ps_full <- cbind(ps, pJ)
  
  if (any(ps_full <= 0)) {
    stop("Some probabilities are estimated to be zero")
  }
  
  loglik <- sum(Y*log(ps_full))
  
  
  return(-loglik)
}

