# adapted from github.com/statdivlab/radEmu, included in raoBust to avoid dependence on radEmu
# helper function for multinomial logistic regression with Firth penalty
get_augmentations <- function(G, Y, B, z) {
  
  n <- length(z)
  p <- nrow(B)
  J <- ncol(Y)
  
  # long version of beta 
  long_beta <- as.vector(B)
  to_remove <- p * (J - 1) + 1:p
  long_beta <- long_beta[-to_remove]
  
  # add z's in to parameter vector 
  theta <- c(long_beta, z)
  
  # expected values 
  log_means <- G %*% theta
  
  # get W and info mat 
  W <- Matrix::Diagonal(x = exp(log_means@x))
  info <- Matrix::crossprod(G, W) %*% G
  info <- Matrix::forceSymmetric(info)
  
  # get symmetric square root of inverse of info matrix
  info_chol <- Matrix::chol(info, pivot = FALSE)
  info_half_inv <- Matrix::solve(info_chol)
  
  # get augmentations
  augs <- matrix(0, nrow = n, ncol = J) 
  for (i in 1:n) {
    G_i <- G[1:J + (i - 1)*J, , drop = FALSE]
    GI_i <- G_i %*% info_half_inv
    augs[i, ] <- Matrix::rowSums(Matrix::Diagonal(x = exp(log_means[(i - 1) * J + 1:J]))
                                 %*% GI_i^2) / 2
  }
  
  return(augs)
  
}