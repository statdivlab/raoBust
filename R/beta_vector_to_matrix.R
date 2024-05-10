#' Create \eqn{\beta} matrix from vector of \eqn{\beta} for \eqn{\beta_k : k \neq j}
#'
#'
#' @param values A vector containing the values for all \eqn{\beta_k}, with \eqn{k \neq j}, as well as all \eqn{\beta_{k0}, for k = 1, \dots, J}.
#' In particular, this vector should be so that the first \eqn{(J-2)(p+1)} entries are \eqn{\beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots, \beta_{(j-1)0}, \beta_{j-1}^{\top},  \beta_{(j+1)0}, \beta_{j+1}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}}.
#' Then, the \eqn{(J-2)(p+1) + 1} entry should be \eqn{\beta_{j0}}.
#' @param p This should be the number of covariates.
#' @param J This should be the number of categories.
#' @param null_j This specifies for which category you have set \eqn{\beta_j = 0}.
#' @param beta_j_null This is the null hypothesized value for \eqn{\beta_j}, which is by default set to be \eqn{\beta_j = 0}.
#'
#' @return The full \eqn{(p+1) \times (J-1)} matrix of \eqn{\beta}
#'
#' @author Shirley Mathur
#'
beta_vector_to_matrix <- function(values, p, J, null_j, beta_j_null = NULL) {
  
  if(is.null(beta_j_null)) {
    beta_j_null <- rep(0,p) #set default value of beta_j_null to 0 if not pre-specified by user
  }
  
  beta_matrix <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
  for (k in 1:(J-1)) {
    if (k < null_j) {
      beta_matrix[ ,k] <- values[((k-1)*(p+1)+1):(k*(p+1))] #populate column k where k < j
    } else if (k == null_j) {
      beta_matrix[,k] <- c(values[(p+1)*(J-2)+1], rep(0,p)) #populate j-th column based on \beta_{j0} value and beta_j_null value
    } else if (k > null_j) {
      beta_matrix[ ,k] <- values[((k-2)*(p+1)+1):((k-1)*(p+1))] #populate column k where k > j (have to shift further down by 1 in index since j not in vector)
    }
  }
  
  return(beta_matrix)
  
}

