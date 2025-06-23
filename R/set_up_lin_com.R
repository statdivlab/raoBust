#' Create matrix to be populated with coefficients for user-specified number of hypotheses.
#'
#' @param J Number of categories.
#'
#' @param p Number of coefficients (excluding intercept).
#' 
#' @param n_hypotheses Number of hypotheses to test.
#'
#' @return The multinomial probabilities for a given value of \eqn{\beta}.
#'
#' @author Shirley Mathur
#'
#'

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
