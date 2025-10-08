#' Compute V matrix contribution to robust score statistic for glm robust score tests.
#'
#' @param indices Indices of observations for cluster.
#' @param model_fits The fitted glm under the null hypothesis for the given indices.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param m Number of parameters fixed under null hypothesis.
#' @param corr_mat Working correlation matrix for model fit.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

V_matrix_contribution <- function(indices, model_fits, yy, xx, m = 1, corr_mat, family, link) {
  V_i <- matrix(NA, nrow = length(model_fits[indices]), ncol = length(model_fits[indices]))
  n_i <- length(indices)
  
  if (family == "poisson" & link == "log") {
    if (n_i > 1) {
      V_i <- diag(sqrt(model_fits[indices])) %*% corr_mat %*% diag(sqrt(model_fits[indices]))
    } else {
      V_i <- model_fits[indices] * corr_mat
    }
  }
  
  if (family == "binomial" & link == "logit") {
    if (n_i > 1) {
      V_i <- diag(sqrt(model_fits[indices]*(1 - model_fits[indices]))) %*% corr_mat %*% diag(sqrt(model_fits[indices]*(1 - model_fits[indices])))
    } else {
      V_i <- model_fits[indices] * corr_mat
    }
  }
  
  if (family == "gaussian") {
    n <- length(yy)
    mod_df <- ncol(xx) - m
    sigma2_tilde <- sum((yy - model_fits)^2)/(n - mod_df)
    if (n_i > 1) {
      V_i <- sqrt(diag(sigma2_tilde, n_i)) %*% corr_mat %*% sqrt(diag(sigma2_tilde, n_i))
    } else {
      V_i <- sigma2_tilde * corr_mat
    }
  }
  
  return (V_i)
}
