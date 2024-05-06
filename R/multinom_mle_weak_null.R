#' Objective function to find Multinomial MLE under weak null ($\beta_j = 0$ for specific j)
#'
#' @param betanonj A vector containing the initial values for all \beta_k, with k \neq j, as well as all \beta_{k0}, for k = 1, \dots, J-1.
#' In particular, this vector should be so that the first (J-2)*(p+1) entries are \beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots, \beta_{(j-1)0}, \beta_{j-1}^{\top},  \beta_{(j+1)0}, \beta_{j+1}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}.
#' @param betaj This is the null hypothesized value for \beta_j, which is by default set to be \beta_j = 0.
#' @param Y This should be the n x J data matrix of outcomes.
#' @param X This should be the n x p design matrix of covariates.
#' @param j This specifies for which category you want to compute the MLE under the constraint that \beta_j = 0.
#'
#' @author Shirley Mathur
#'
#' @return A vector containing the optimal `betanonj` values to maximize the log-likelihood under the null constraint that \beta_j = 0. The components are listed out in the same manner as in the betanonj parameter.
#'
multinom_mle_weak_null <- function(betanonj, betaj = rep(0, p), Y, X, j) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
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


