#' Robust score (Rao) test for Poisson regression
#'
#' @param glm_object The fitted glm under the alternative.
#' @param call_to_model The call used to fit the model. Used internally.
#' @param param which parameter do you want to test?
#' @param id observations with the same id are in the same cluster
#'
#' @importFrom stats glm poisson model.matrix glm.fit pchisq
#'
#' @export
robust_score_test <- function(glm_object, call_to_model, param = 1,
                              id = NA) {

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

  if (is.factor(id)) {

    ## assume exchangeable

    ids <- unique(id)
    Umatrices <- list()
    Amatrices <- list()


    for (ii in 1:length(ids)) {
      indices <- which(id == ids[ii])
      n_i <- length(indices)
      xxi <- xx[indices, , drop = FALSE]
      model0_fits_i <- model0_fits[indices]

      Di <- matrix(NA, nrow = pp, ncol = n_i)
      for (j in 1:n_i) {
        for (k in 1:pp) {
          Di[k, j] <- xxi[j, k] * model0_fits_i[j]
        }
      }

      corr_matrix <- matrix(rep(glm_object$geese$alpha, n_i^2), nrow = n_i)
      diag(corr_matrix) <- 1
      if (n_i > 1) {
        Vi <- diag(sqrt(model0_fits_i)) %*% corr_matrix %*% diag(sqrt(model0_fits_i))
      } else {
        Vi <- sqrt(model0_fits_i) * corr_matrix * sqrt(model0_fits_i)
      }

      Si <- yy[indices] - model0_fits_i

      Umatrices[[ii]] <- Di %*% solve(Vi) %*% matrix(Si, ncol = 1)
      Amatrices[[ii]] <- Di %*% solve(Vi) %*% t(Di)

    }

    u_tilde_sum <- Reduce("+", Umatrices) # p x 1

    aa0 <- Reduce("+", Amatrices) ### p x p
    aa0_11 <- aa0[setdiff(1:pp, param), setdiff(1:pp, param)]
    aa0_22 <- aa0[param, param]
    aa0_21 <- aa0[param, setdiff(1:pp, param)]
    c_tilde <- cbind(-aa0_21 %*% solve(aa0_11), diag(rep(1, pp0))) # 1 x p

    test_stat <- c(t(c_tilde %*% u_tilde_sum) %*% ## (1 x p) x (p x 1)
                     solve(c_tilde %*% (Reduce("+", lapply(FUN = function(x) { x %*% t(x) }, Umatrices))) %*% t(c_tilde)) %*% # (1 x p) x ((p x n) x (n x p)) x (p x 1)
                     (c_tilde %*% u_tilde_sum))

  } else if (is.na(id)) {

    score_contribution <- function(i, model_fits) { ### output is p x 1
      matrix(model_fits[i]^2 * (yy[i] - model_fits[i]) * xx[i, ], ncol = 1)
    }

    fisher_info_contribution <- function(i, model_fits) { ### output is p x p... Di^T Vi Di
      model_fits[i] * xx[i, ] %*% t(xx[i, ])
    }

    u_tilde <- sapply(1:nn, score_contribution, model_fits = model0_fits) ### p x n
    u_tilde <- rbind(u_tilde[-param,], u_tilde[param,]) ### reorder U_tilde the way A is reordered

    aa0 <- Reduce("+", sapply(1:nn, fisher_info_contribution, model_fits = model0_fits, simplify=F)) ### p x p
    aa0_11 <- aa0[setdiff(1:pp, param), setdiff(1:pp, param)]
    aa0_22 <- aa0[param, param]
    aa0_21 <- aa0[param, setdiff(1:pp, param)]
    c_tilde <- cbind(-aa0_21 %*% solve(aa0_11), diag(rep(1, pp0))) # 1 x p

    test_stat <- c(t(c_tilde %*% matrix(rowSums(u_tilde), ncol = 1)) %*% ## (1 x p) x (p x 1)
                     solve(c_tilde %*% (u_tilde %*% t(u_tilde)) %*% t(c_tilde)) %*% # (1 x p) x ((p x n) x (n x p)) x (p x 1)
                     (c_tilde %*% matrix(rowSums(u_tilde), ncol = 1)))
  } else {

    stop("unsure what correlation structure is")

  }

  1 - pchisq(test_stat, df = pp0)

}
