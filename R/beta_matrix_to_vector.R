beta_matrix_to_vector <- function(beta_vector, null_j = NULL) {
  
  # if (is.null(null_j)) {
  #   null_j <- rep(0, p)
  # }
  #
  # stopifnot(length(values) == ((p + 1) * (J - 1)) - p)
  #
  # betas_no_jth_col <- matrix(values[-null_j], nrow = p + 1, ncol = J - 2, byrow = TRUE)
  #
  # cbind(betas_no_jth_col[, 1:(null_j - 1)],
  #       c(values[null_j], rep(0, p)),
  #       betas_no_jth_col[, null_j: (J-2)])
  stop("this function needs to be implemented")
}
