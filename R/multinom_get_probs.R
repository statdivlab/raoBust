#' Compute multinomial probabilities for a given value of model parameters.
#'
#' @param Y This should be the \eqn{n \times J} data matrix of outcomes.
#' @param X This should be the \eqn{n \times p} design matrix of covariates.
#' @param beta This should be the \eqn{(p+1) \times (J-1) \beta} matrix.
#'
#' @return The multinomial probabilities for a given value of \eqn{\beta}.
#'
#' @author Shirley Mathur
#'
#'
multinom_get_probs <- function(X, Y, beta) {
  #terms necessary for computation of probabilities
  Xaug <- cbind(1, X)
  XaugBeta <- Xaug %*% beta
  pJ <- (1 + rowSums(exp(XaugBeta)))^(-1)
  ps <- as.vector(pJ)*exp(XaugBeta)
  ps_full <- cbind(ps, pJ)
  
  return(ps_full)
}
