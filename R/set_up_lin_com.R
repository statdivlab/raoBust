#' Create matrix to be populated with coefficients for user-specified number of hypotheses.
#'
#' @param J Number of categories.
#'
#' @param p Number of coefficients (excluding intercept).
#' 
#' @param n_hypotheses Number of hypotheses to test.
#'
#' @return An 'A' matrix with `n_hypotheses` rows and `(J-1) x (p+1)` columns with the various coefficient and category combinations as column names.
#' The matrix is filled with all 0's, but can be modified by the user to reflect the linear combinations they wish to test.
#'
#' @author Shirley Mathur
#'
#'@export

set_up_lin_com <- function (J, p, n_hypotheses) {
  
  num_params <- (p+1)*(J-1)
  A <- matrix(0, ncol = num_params, nrow = n_hypotheses)
  colnames(A) <- paste(rep("k", num_params),
                       rep(1:(p+1) - 1, J-1),
                       rep("j", num_params),
                       rep(1:(J-1), each = p+1),
                       sep = "_")
  
  return(A)
  
}
