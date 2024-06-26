#' Compute Information matrix for model parameters.
#'
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param probs This should be the \eqn{n \times J} matrix of estimated probabilities for each sample.
#'
#' @return The Fisher information matrix for model parameters.
#'
#' @author Shirley Mathur
#'
#'
multinom_info_mat <- function(X, Y, probs) {
  #get relevant quantities necessary to compute score vector
  p <- ncol(X)
  J <- ncol(Y)
  n <- nrow(Y)
  N <- matrix(rowSums(Y), nrow = n) #totals by sample
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  
  info_mat <- matrix(data = rep(0, (p+1)^2*(J-1)^2), ncol = (p+1)*(J-1)) #initialize matrix for information matrix
  for (k in 1:(J-1)) {
    for (l in 1:k) {
      if (l != k) {
        info_mat[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <-  t(Xaug) %*% (as.vector(-probs[ ,k]*probs[ ,l]*N)*Xaug)
        info_mat[((k-1)*(p+1) + 1):(k*(p+1)), ((l-1)*(p+1) + 1):(l*(p+1))] <- t(info_mat[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))])
      } else if (l == k) {
        info_mat[((k-1)*(p+1) + 1):(k*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (as.vector(-probs[ ,k]*(-1 + probs[ ,k])*N)*Xaug)
      }
    }
  }
  
  
  return(info_mat)
}
