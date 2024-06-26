#' Compute score for model parameters.
#'
#' @param Y This should be the \eqn{n \times J} data matrix of outcomes.
#' @param X This should be the \eqn{n \times p} design matrix of covariates.
#' @param probs This should be the \eqn{n \times J} matrix of estimated probabilities for each sample.
#'
#' @return The score vector for the model parameters.
#'
#' @author Shirley Mathur
#'
#'
multinom_score_vector <- function(X, Y, probs) {
  #get relevant quantities necessary to compute score vector
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)
  N <- matrix(rowSums(Y), nrow = n) #totals by sample
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  
  score <- matrix(data = rep(0, (p+1)*(J-1)), ncol = 1) #initialize matrix for score
  for (k in 1:(J-1)) {
    score[((k-1)*(p+1) + 1):(k*(p+1)) ,] <- colSums(as.vector(Y[ ,k] - probs[ ,k]*N)*Xaug)
  }
  
  return(score)
}
