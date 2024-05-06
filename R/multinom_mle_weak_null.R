#' Objective function to find Multinomial MLE under weak null (\eqn{\beta_j = 0} for specific j)
#'
#'
#' @param betanonj A vector containing the initial values for all \eqn{\beta_k}, with \eqn{k \neq j}, as well as all \eqn{\beta_{k0}, for k = 1, \dots, J-1}.
#' In particular, this vector should be so that the first \eqn{(J-2)(p+1)} entries are \eqn{\beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots, \beta_{(j-1)0}, \beta_{j-1}^{\top},  \beta_{(j+1)0}, \beta_{j+1}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}}.
#' @param betaj This is the null hypothesized value for \eqn{\beta_j}, which is by default set to be \eqn{\beta_j = 0}.
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param j This specifies for which category you want to compute the MLE under the constraint that \eqn{\beta_j = 0}.
#'
#' @return A vector containing the optimal betanonj values to maximize the log-likelihood under the null constraint that \eqn{\beta_j = 0}. The components are listed out in the same manner as in the betanonj parameter.
#'
#'
#' @author Shirley Mathur
#'
#'
multinom_mle_weak_null <- function(betanonj, betaj = NULL, Y, X, j) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  p <- ncol(X)
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  if (is.null(betaj)) {
    betaj <- rep(0, p)
  }
  beta <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
  for (k in 1:(J-1)) {
    if (k < j) {
      beta[ ,k] <- betanonj[((k-1)*(p+1)+1):(k*(p+1))]
    } else if (k == j) {
      beta[,k] <- c(betanonj[(p+1)*(J-2)+1], betaj)
    } else if (k > j) {
      beta[ ,k] <- betanonj[((k-2)*(p+1)+1):((k-1)*(p+1))]
    }
  }

  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  XaugBeta <- Xaug %*% beta
  pJ <- (1 + rowSums(exp(XaugBeta)))^(-1)
  ps <- as.vector(pJ)*exp(XaugBeta)
  ps_full <- cbind(ps, pJ)

  loglik <- sum(Y*log(ps_full))

  return(-loglik)
}


