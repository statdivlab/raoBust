#' Objective function to find Multinomial MLE under global null ($\beta_1 = beta_2 = \dots = \beta_{J-1} = 0$)
#'
#' @param betanots A vector containing the initial values for all \beta_{k0}, for k = 1, \dots, J-1.
#' In particular, this vector should be so that the entries are ordered as \beta_{10}, \beta_{20}, \dots, \beta_{(J-1)0}.
#' @param Y This should be the n x J data matrix of outcomes.
#' @param X This should be the n x p design matrix of covariates.
#' @author Shirley Mathur
#' @return A vector containing the optimal values for `betanots`  to maximize the log-likelihood under the null constraint that \beta_1 = beta_2 = \dots = \beta_{J-1} = 0.
#' The components are listed out in the same manner as in the `betanots` parameter.
multinom_mle_global_null <- function(betanots, Y, X) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample

  pJ <- (1 + sum(exp(betanots)))^(-1)
  ps <- pJ*exp(betanots)
  ps_full <- matrix(rep(c(ps, pJ), n), nrow = n, byrow = TRUE)

  loglik <- sum(Y*log(ps_full))


  return(-loglik)
}

