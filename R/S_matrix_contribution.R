#' Compute S matrix contribution for a given cluster to robust score statistic for glm robust score tests.
#'
#' @param indices Indices of observations in cluster.
#' @param model_fits The fitted glm under the null hypothesis.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

S_matrix_contribution <- function(indices, model_fits, yy, xx, family, link) {

  S_i <- matrix(yy[indices] - model_fits[indices], ncol = 1)

  return (S_i)
}
