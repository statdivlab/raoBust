#' Robust score (Rao) test for Poisson regression
#'
#' @param glm_object The fitted glm under the alternative.
#' @param call_to_model The call used to fit the model. Used internally.
#' @param param which parameter do you want to test?
#'
#' @importFrom stats glm poisson model.matrix glm.fit pchisq
#'
#' @export
robust_score_test <- function(glm_object, call_to_model, param = 1) {

  # stop("no")
  model1 <- glm_object
  model1_estimates <- model1$coef
  model1_fits <- model1$fitted.values

  xx <- model.matrix(model1)
  yy <- model1$y
  nn <- nrow(xx)
  pp <- ncol(xx)
  xx0 <- xx[ , -param]
  pp0 <- length(param)

  # stop("no")
  withCallingHandlers({
    model0 <- glm.fit(x = xx0,
                      y = yy,
                      intercept = FALSE,
                      offset = with(model1$data, eval(call_to_model$offset)),
                      weights = with(model1$data, eval(call_to_model$weights)),
                      family = eval(call_to_model$family))

  }, warning = function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })

  model0_fits <- model0$fitted.values
  fitted_cols0 <- matrix(rep(model0_fits, each = pp - 1), nrow = nn, ncol = pp - 1, byrow = T)

  score_contribution <- function(i, model_fits) {
    matrix(model_fits[i]^2 * (yy[i] - model_fits[i]) * xx[i, ], ncol = 1)
  }

  u_tilde <- sapply(1:nn, score_contribution, model_fits = model0_fits)

  fisher_info_contribution <- function(i, model_fits) {
    model_fits[i] * xx[i, ] %*% t(xx[i, ])
  }

  aa0 <- Reduce("+", sapply(1:nn, fisher_info_contribution, model_fits = model0_fits, simplify=F))
  aa0_11 <- aa0[setdiff(1:pp, param), setdiff(1:pp, param)]
  aa0_22 <- aa0[param, param]
  aa0_21 <- aa0[param, setdiff(1:pp, param)]
  c_tilde <- cbind(-aa0_21 %*% solve(aa0_11), diag(rep(1, pp0))) ## maybe?

  test_stat <- c(t(c_tilde %*% matrix(rowSums(u_tilde), ncol = 1)) %*%
                   solve(c_tilde %*% (u_tilde %*% t(u_tilde)) %*% t(c_tilde)) %*%
                   (c_tilde %*% matrix(rowSums(u_tilde), ncol = 1)))

  1 - pchisq(test_stat, df = pp0)

}
