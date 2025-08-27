#' Compute D matrix contribution for a given cluster to robust score statistic for glm robust score tests.
#'
#' @param indices Indices of observations in cluster.
#' @param model_fits The fitted glm under the null hypothesis.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

D_matrix_contribution <- function(indices, model_fits, yy, xx, family, link) {
  D_i <- matrix(NA, nrow = ncol(xx), ncol = length(model_fits))
  
  if (family == "poisson" & link == "log") {
    D_i <- t(model_fits[indices]*xx[indices,])
  }
  
  if (family == "binomial" & link == "logit") {
    D_i <- t(model_fits[indices]*(1-model_fits[indices])*xx[indices,])
  }
  
  if (family == "gaussian" & link == "identity") {
    D_i <- t(xx[indices,])
  }
  
  return (D_i)
}
