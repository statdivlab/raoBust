#' Robust score (Rao) tests with finite-sample correction
#'
#' The default behavior is to do no finite sample correction (for the covariance of the score) for correlated data,
#' but to do it for uncorrelated data. This choice performed best under a small simulation study. See Guo et al
#' for details; the proposed modification is in equation 20.
#'
#' @param glm_object The fitted glm under the alternative.
#' @param call_to_model The call used to fit the model. Used internally.
#' @param param which parameter do you want to test? Used internally.
#' @param id observations with the same id are in the same cluster
#'
#' @references
#' Guo, X., Pan, W., Connett, J. E., Hannan, P. J., & French, S. A. (2005).
#' Small-sample performance of the robust score test and its modifications in
#' generalized estimating equations. *Statistics in Medicine, 24*(22), 3479–3495.
#' Wiley Online Library. <doi:10.1002/sim.2161>

#'
#' @importFrom stats glm poisson model.matrix glm.fit pchisq
#'
#' @export
robust_score_test <- function(glm_object, call_to_model, param = 1,
                              id = NA) {

  model1 <- glm_object
  model1family <- glm_object$family$family
  model1link <- glm_object$family$link
  model1_estimates <- model1$coef
  model1_fits <- model1$fitted.values

  xx <- model.matrix(model1)
  yy <- model1$y
  nn <- nrow(xx)
  pp <- ncol(xx)
  xx0 <- xx[ , -param, drop = FALSE]
  pp0 <- length(param)



  # stop("no")
  withCallingHandlers({
    if (model1family == "gaussian" & model1link == "log") {
      model0 <- glm.fit(x = xx0,
                        y = yy,
                        intercept = FALSE,
                        offset = with(model1$data, eval(call_to_model$offset)),
                        weights = with(model1$data, eval(call_to_model$weights)),
                        family = eval(call_to_model$family),
                        start = rep(1, ncol(xx0)))
    } else {
      model0 <- glm.fit(x = xx0,
                        y = yy,
                        intercept = FALSE,
                        offset = with(model1$data, eval(call_to_model$offset)),
                        weights = with(model1$data, eval(call_to_model$weights)),
                        family = eval(call_to_model$family))
    }
    

  }, warning = function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })
  model0_fits <- model0$fitted.values

  if (length(id) == nn) {

    ## assume exchangeable

    ids <- unique(id)
    Umatrices <- list()
    Amatrices <- list()


    for (ii in 1:length(ids)) {
      indices <- which(id == ids[ii])
      n_i <- length(indices)
      xxi <- xx[indices, , drop = FALSE]
      model0_fits_i <- model0_fits[indices]

      Di <- D_matrix_contribution(indices = indices, model_fits = model0_fits,
                                  yy = yy, xx = xx, family = model1family, link = model1link)

      corr_matrix <- matrix(rep(glm_object$geese$alpha, n_i^2), nrow = n_i)
      diag(corr_matrix) <- 1
      Vi <- V_matrix_contribution(indices = indices, model_fits = model0_fits, corr_mat = corr_matrix,
                                  yy = yy, xx = xx, family = model1family, link = model1link)

      Si <- S_matrix_contribution(indices = indices, model_fits = model0_fits,
                                  yy = yy, xx = xx, family = model1family, link = model1link)

      ## Recall solve(x1, x2) is the same as solve(x1) %*% x2
      Umatrices[[ii]] <- Di %*% solve(Vi, Si)
      Amatrices[[ii]] <- Di %*% solve(Vi, t(Di))

    }

    ### Order of parameters matters. Set up permutation matrix
    perm_idx <- c(setdiff(seq_len(pp), param), param)
    Pi <- diag(pp)[perm_idx, , drop = FALSE]  # p x p permutation

    ### Compute U and reorder
    u_tilde_sum <- Pi %*% Reduce("+", Umatrices) # p x 1
    # u_tilde_sum <- matrix(c(u_tilde_sum[-param,1], u_tilde_sum[param,1]), ncol = 1)

    ### Compute B and reorde
    # cov_U_hat <- (Reduce("+", lapply(FUN = function(x) { x %*% t(x) }, Umatrices)))
    cov_U_hat <- Reduce("+", lapply(Umatrices, function(u) {
      up <- Pi %*% u
      up %*% t(up)
    }))

    ## finite sample correction from Guo et al, addresses conservatism in small samples and improves power
    ## under simulation with random effects, the test was actually anticonservative for small samples
    ## for this reason, we don't do this by default
    # if (length(ids) > 1L) {
    #   cov_U_hat <- (length(ids) - 1) / length(ids) * cov_U_hat
    # }

    ### Compute A
    aa0 <- Reduce("+", Amatrices) ### p x p
    ## Reorder A
    aa0_11 <- aa0[setdiff(seq_len(pp), param), setdiff(seq_len(pp), param), drop = FALSE]
    aa0_21 <- aa0[param, setdiff(seq_len(pp), param), drop = FALSE]
    # c_tilde <- cbind(-aa0_21 %*% solve(aa0_11), diag(rep(1, pp0))) # 1 x p
    c_tilde <- cbind(-t(solve(aa0_11, t(aa0_21))), diag(1, pp0))


    ## Reorder B
    test_stat <- c(t(c_tilde %*% u_tilde_sum) %*% ## (1 x p) x (p x 1)
                     solve(c_tilde %*% cov_U_hat %*% t(c_tilde), # (1 x p) x ((p x n) x (n x p)) x (p x 1)
                           (c_tilde %*% u_tilde_sum)))

  } else if (length(id) != nn) {

    u_tilde <- sapply(1:nn, score_contribution, model_fits = model0_fits,
                      family = model1family, link = model1link, xx = xx, yy = yy) ### p x n

    ## Reorder u tilde to match a
    u_tilde <- rbind(u_tilde[-param,], u_tilde[param,])
  
    aa0 <- Reduce("+", sapply(1:nn, fisher_info_contribution, simplify=F,  model_fits = model0_fits,
                              family = model1family, link = model1link, xx = xx, yy = yy, m = pp0)) ### p x p
    
    aa0_11 <- aa0[setdiff(seq_len(pp), param), setdiff(seq_len(pp), param), drop = FALSE]
    aa0_21 <- aa0[param, setdiff(seq_len(pp), param), drop = FALSE]
    # c_tilde <- cbind(-aa0_21 %*% solve(aa0_11), diag(rep(1, pp0))) # 1 x p
    c_tilde <- cbind(-t(solve(aa0_11, t(aa0_21))), diag(1, pp0))

    ## (nn - 1) / nn is finite sample correction from Guo et al, addresses conservatism in small samples and improves power
    test_stat <- c(t(c_tilde %*% matrix(rowSums(u_tilde), ncol = 1)) %*% ## (1 x p) x (p x 1)
                     solve(c_tilde %*% ((nn - 1) / nn * u_tilde %*% t(u_tilde)) %*% t(c_tilde), # (1 x p) x ((p x n) x (n x p)) x (p x 1)
                           (c_tilde %*% matrix(rowSums(u_tilde), ncol = 1))))

  } else {

    stop("unsure what correlation structure is")

  }

  1 - pchisq(test_stat, df = pp0)

}
