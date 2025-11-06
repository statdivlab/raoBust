#' Robust score (Rao) tests for multinomial regression.
#'
#' @param X A \eqn{n x p} design matrix of covariates.
#' @param Y A \eqn{n x J} data matrix of outcomes.
#' @param formula a one-sided formula specifying the form of the mean model to be fit (use with \code{data} argument if \code{X} is not included)
#' @param data a dataframe with \eqn{n} rows containing variables given in \code{formula} (use with \code{formula} argument if \code{X} is not included)
#' @param strong If FALSE, this function will compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0} for the length \eqn{p} vector \eqn{\beta_j}.
#' If TRUE, this function instead computes the robust score statistic to test the strong null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0} for all length \eqn{p} vectors \eqn{\beta_j}, \eqn{j\in\{1,\ldots,J-1\}}.
#' Default is FALSE.
#' @param j If `strong` is FALSE, this argument must be supplied. This gives the category \eqn{j} in the weak null hypothesis that \eqn{\beta_j = 0}.
#' @param all_score If TRUE, score tests for each individual covariate and category pair (i.e. null that \eqn{\beta_{jk} = 0} for each category \eqn{j = 1, \dots, J-1} and each covariate \eqn{k = 1, \dots, p} pair) will be
#' run and reported in output coefficient table. Default is FALSE.
#' @param penalty If TRUE will apply a Firth penalty to estimation under the alternative and under the null. Defaults to FALSE (ask Amy her preference)
#' @param pseudo_inv Use pseudo inverse for inverted portion of the robust score test to avoid issues with nearly singular matrices.
#'
#' @return The robust score test statistic for the specified hypothesis test. A list including the test statistic, p-value,
#' estimated parameters under the null hypothesis, and estimated parameters under the alternative hypothesis.
#'
#' @importFrom stats nlm optim formula qnorm
#'
#' @author Shirley Mathur
#'
#' @export
multinom_test <- function(X = NULL, Y, formula = NULL, data = NULL,
                          strong = FALSE, j = NULL, all_score = FALSE, penalty = FALSE, pseudo_inv = FALSE) {

  #record function call
  cl <- match.call()

  # if X is null and formula and data are provided, get design matrix
  if (is.null(X)) {
    if (is.null(formula) | is.null(data)) {
      stop("If design matrix X not provided, both formula and data containing
covariates in formula must be provided.")
    }
    X <- model.matrix(formula, data)
  }

  # if X has intercept column (or a column of constant values across all observations), remove it
  constant_cols <- apply(X, 2, function(x) {all(x == mean(x, na.rm = TRUE))})
  if (sum(constant_cols) > 0) {
    X <- X[, !constant_cols, drop = FALSE]
  }

  # check that X and Y have the same number of rows, if not throw error
  if (nrow(X) != nrow(Y)) {
    stop("Please make sure that X and Y have the same number of observations.")
  }

  # check that if X and Y have rownames they are the same, if not throw warning
  if (!is.null(rownames(X)) & !is.null(rownames(Y))) {
    if (!(all.equal(rownames(X), rownames(Y)) == TRUE)) {
      warn("X and Y both have rownames, and they do not match. Please check that the rows of X and Y correspond to the
         same observations.")
    }
  }

  #get n, p, J values (used throughout rest of the function to compute relevant quantities)
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)

  #compute other necessary quantities needed for computations
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix

  #compute test statistics under weak null that \beta_j = 0 for user-specified j
  if (strong == FALSE) {
    #report error if user specified marginal test but does not supply an argument for j
    if(is.null(j)) {
      stop("If testing under the weak null hypothesis with `strong = FALSE`, you must include the argument `j`.")
    }

    #get beta mle under weak null constraint, first setting up initial value of beta and then optimizing
    beta_init <- matrix(1, nrow = p + 1, ncol = J-1)
    beta_init[(2:(p+1)), j] <- 0
    if (!penalty) {
      beta_null1mle <- multinom_fisher_scoring(beta_init, X, Y, strong = FALSE, null_j = j, pseudo_inv = pseudo_inv)
    } else {
      beta_null1mle <- multinom_penalized_estimation(beta_init, X, Y, strong = FALSE, null_j = j, pseudo_inv = pseudo_inv)
    }

    #jacobian of function of parameter h(\beta) = 0
    H1 <- matrix(data = rep(0, p*(p+1)*(J-1)), ncol = (p+1)*(J-1), nrow = p)
    H1[1:p, ((j-1)*(p+1) + 2):(j*(p+1))] <- diag(nrow = p, ncol = p)

    #terms necessary for computation of S, I, D matrices
    ps_full <- multinom_get_probs(X, Y, beta_null1mle)


    #score of \beta evaluated at mle under null constraint
    S1 <- multinom_score_vector(X, Y, ps_full)

    #information matrix evaluated at mle under null constraint
    I1 <- multinom_info_mat(X, Y, ps_full)

    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D1 <- matrix(data = rep(0, length(S1)^2), ncol = length(S1))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S1_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S1_l <- as.vector(Y[ ,l] - ps_full[ ,l]*N)*Xaug
        D1[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S1_k)%*%S1_l
      }
    }

    the_mle <- beta_null1mle
    the_df <- p


    #compute statistic!
    if (!pseudo_inv) {
      ## Recall solve(x1, x2) is the same as solve(x1) %*% x2
      T_GS<- tryCatch({as.numeric(t(S1) %*% solve(I1, t(H1)) %*%
                                    solve(H1 %*% solve(I1, D1) %*% solve(I1, t(H1)),
                                          H1 %*% solve(I1, S1)))},
                      error = function(cond) {
                        print(cond)
                        return(NA)
                      })
    } else {
      T_GS<- tryCatch({as.numeric(t(S1) %*% MASS::ginv(I1) %*% t(H1) %*%
                                    (MASS::ginv(H1 %*% MASS::ginv(I1) %*% D1 %*% MASS::ginv(I1) %*% t(H1))) %*%
                                    H1 %*% MASS::ginv(I1) %*% S1)},
                      error = function(cond) {
                        print(cond)
                        return(NA)
                      })
    }


  } else {
    #initialize value for beta_k0 for all k = 1, \dots, J-1
    #betanots <- rep(1, J-1)

    #jacobian of function of parameter h(\beta) = 0
    H2 <- matrix(data = rep(0, p*(J-1)*(p+1)*(J-1)), ncol = (p+1)*(J-1), nrow = p*(J-1))
    for (k in 1:(J-1)) {
      H2[ ((k-1)*p + 1):(k*p), ((k-1)*(p+1) + 2):(k*(p+1))] <- diag(nrow = p, ncol = p)
    }

    #compute mle under null constraint
    beta_init <- matrix(0, nrow = p + 1, ncol = J-1)
    beta_init[1,] <- 1
    if (!penalty) {
      beta_null2mle <- multinom_fisher_scoring(beta_init, X, Y, strong = TRUE, pseudo_inv = pseudo_inv)
    } else {
      beta_null2mle <- multinom_penalized_estimation(beta_init, X, Y, strong = TRUE, pseudo_inv = pseudo_inv)
    }

    ## AW TODO below not needed?
    # beta_null2mle <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
    # beta_null2mle[1, ] <- betanots_null2mle


    #terms necessary for computation of S, I, D matrices
    ps_full <- multinom_get_probs(X, Y, beta_null2mle)


    #score of \beta evaluated at mle under null constraint
    S2 <- multinom_score_vector(X, Y, ps_full)

    #information matrix evaluated at mle under null constraint
    I2 <- multinom_info_mat(X, Y, ps_full)

    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D2 <- matrix(data = rep(0, length(S2)^2), ncol = length(S2))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S2_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S2_l <- as.vector(Y[ ,l] - ps_full[ ,l]*N)*Xaug
        D2[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S2_k)%*%S2_l
      }
    }

    the_mle <- beta_null2mle
    the_df <- p*(J-1)

    #compute statistic!
    if (!pseudo_inv) {
      T_GS <- tryCatch({as.numeric(t(S2) %*% solve(I2, t(H2)) %*%
                                     solve(H2 %*% solve(I2, D2) %*% solve(I2, t(H2)),
                                           H2 %*% solve(I2, S2)))},
                       error = function(cond) {
                         warning("Test statistic cannot be calculated due to the error printed above.")
                         print(cond)
                         return(NA)
                       })
    } else {
      T_GS <- tryCatch({as.numeric(t(S2) %*% MASS::ginv(I2) %*% t(H2) %*%
                                     (MASS::ginv(H2 %*% MASS::ginv(I2) %*% D2 %*% MASS::ginv(I2) %*% t(H2))) %*%
                                     H2 %*% MASS::ginv(I2) %*% S2)},
                       error = function(cond) {
                         warning("Test statistic cannot be calculated due to the error printed above.")
                         print(cond)
                         return(NA)
                       })
    }



  }


  #run score test for each coefficient (all covariate and category pairs) if desired
  score_test_pvals <- matrix(NA, nrow = p+1, ncol = J-1)
  if (all_score) {
    for (j in 1:(J-1)) {
      for (m in 1:(p+1)) {
        #jacobian of function of parameter h(\beta) = 0
        H3 <- matrix(data = rep(0, (p+1)*(J-1)), nrow = 1)
        H3[1,(j-1)*(p+1)+m] <- 1

        #compute mle under null constraint
        beta_init <- matrix(1, nrow = p + 1, ncol = J-1)
        beta_init[m,j] <- 0
        if (!penalty) {
          beta_null3mle <- multinom_fisher_scoring(beta_init, X, Y, j_ind = j, k_ind = m, pseudo_inv = pseudo_inv)
        } else {
          beta_null3mle <- multinom_penalized_estimation(beta_init, X, Y, j_ind = j, k_ind = m, pseudo_inv = pseudo_inv)
        }

        #terms necessary for computation of S, I, D matrices
        ps_full <- multinom_get_probs(X, Y, beta_null3mle)


        #score of \beta evaluated at mle under null constraint
        S3 <- multinom_score_vector(X, Y, ps_full)

        #information matrix evaluated at mle under null constraint
        I3 <- multinom_info_mat(X, Y, ps_full)

        #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
        D3 <- matrix(data = rep(0, length(S3)^2), ncol = length(S3))
        for (k in 1:(J-1)) {
          for (l in 1:(J-1)) {
            S3_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
            S3_l <- as.vector(Y[ ,l] - ps_full[ ,l]*N)*Xaug
            D3[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S3_k)%*%S3_l
          }
        }

        temp_mle <- beta_null3mle

        #compute statistic!
        if (!pseudo_inv) {
          temp_GS <- tryCatch({as.numeric(t(S3) %*% solve(I3, t(H3)) %*%
                                            solve(H3 %*% solve(I3, D3) %*% solve(I3, t(H3)),
                                                  H3 %*% solve(I3, S3)))},
                              error = function(cond) {
                                warning("Test statistic cannot be calculated due to the error printed above.")
                                print(cond)
                                return(NA)
                              })
        } else {
          temp_GS <- tryCatch({as.numeric(t(S3) %*% MASS::ginv(I3) %*% t(H3) %*%
                                            (MASS::ginv(H3 %*% MASS::ginv(I3) %*% D3 %*% MASS::ginv(I3) %*% t(H3))) %*%
                                            H3 %*% MASS::ginv(I3) %*% S3)},
                              error = function(cond) {
                                warning("Test statistic cannot be calculated due to the error printed above.")
                                print(cond)
                                return(NA)
                              })
        }

        #get pval for test and record
        score_test_pvals[m,j] <- pchisq(temp_GS, df = 1, lower.tail = FALSE)
      }
    }

  }


  #compute mle under alternative
  beta_alt <- matrix(-0.02, nrow = p + 1, ncol = J-1)
  if (!penalty) {
    mle_alt <- tryCatch({multinom_fisher_scoring(beta_alt, X, Y, null = FALSE, pseudo_inv = pseudo_inv)},
                        error = function(cond) {
                          warning("MLE under alternative hypothesis cannot be calculated.")
                          print(cond)
                          return(NA)
                        })
  } else {
    mle_alt <- tryCatch({multinom_penalized_estimation(beta_alt, X, Y, null = FALSE, pseudo_inv = pseudo_inv)},
                        error = function(cond) {
                          warning("MLE under alternative hypothesis cannot be calculated.")
                          print(cond)
                          return(NA)
                        })
  }

  #compute robust Wald SE's
  #terms necessary for computation of S, I, D matrices
  ps_alt <-   multinom_get_probs(X, Y, mle_alt)

  #score of \beta evaluated at mle under alternative
  S_alt <- multinom_score_vector(X, Y, ps_alt)

  #information matrix evaluated at mle under alternative
  I_alt <- multinom_info_mat(X, Y, ps_alt)

  #D matrix (sum of S_iS_i^T for all i = 1, \dots, n) for mle under alternative
  D_alt <- matrix(data = rep(0, length(S_alt)^2), ncol = length(S_alt))
  for (k in 1:(J-1)) {
    for (l in 1:(J-1)) {
      S_alt_k <- as.vector(Y[ ,k] - ps_alt[ ,k]*N)*Xaug
      S_alt_l <- as.vector(Y[ ,l] - ps_alt[ ,l]*N)*Xaug
      D_alt[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S_alt_k)%*%S_alt_l
    }
  }

  robust_wald_cov <- solve(I_alt, D_alt) %*% solve(I_alt)
  robust_wald_se <- matrix(sqrt(diag(robust_wald_cov)), nrow = p+1, ncol = J-1, byrow = FALSE)

  #make table of coefficients
  cat_labs <- 1:ncol(Y)
  if(!is.null(colnames(Y))) {
    cat_labs <- colnames(Y)
  }

  coef_tab <- data.frame("Category" = rep(cat_labs[1:J-1], each = p+1),
                         "Covariate" = rep(1:(p+1), J-1),
                         "Estimate" = rep(NA, (p+1)*(J-1)),
                         "Robust Std Error" = rep(NA, (p+1)*(J-1)),
                         "Lower 95% CI" = rep(NA, (p+1)*(J-1)),
                         "Upper 95% CI" = rep(NA, (p+1)*(J-1)),
                         "Robust Wald p" = rep(NA, (p+1)*(J-1)),
                         "Robust Score p" = rep(NA, (p+1)*(J-1)),
                         check.names = FALSE)

  #set covariate terms appropriately if formula was used in call
  if (!is.null(formula)) {

    #get names of variables used for model fit, and then use these to populate covariate column of output table
    coef_names <- c(labels(stats::terms(stats::as.formula(formula), data = data)))
    coef_tab$Covariate <- rep(c("(intercept)",coef_names), J-1)

  } else {
    if (is.null(colnames(X))) {coef_names <- 1:ncol(X)} else {coef_names <- colnames(X)}
    coef_tab$Covariate <- rep(c("(intercept)", coef_names), J-1)
  }

  #populate estimate, se, Wald p, and lower and upper columns of output table with appropriate quantities
  coef_tab$Estimate <- c(mle_alt)
  coef_tab$'Robust Std Error' <- c(robust_wald_se)
  coef_tab$'Robust Wald p' <- pchisq((coef_tab$'Estimate'/coef_tab$'Robust Std Error')^2, df=1, lower.tail = FALSE)
  coef_tab$'Lower 95% CI' <- coef_tab$Estimate + qnorm(0.025)*coef_tab$'Robust Std Error'
  coef_tab$'Upper 95% CI' <- coef_tab$Estimate + qnorm(0.975)*coef_tab$'Robust Std Error'
  coef_tab$'Robust Score p' <- c(score_test_pvals)

  #sort table by magnitude of effect size
  coef_tab <- coef_tab[order(abs(coef_tab$Estimate), decreasing = TRUE), ]

  result <- list("call" = cl,
                 "design_matrix" = Xaug,
                 "test_stat" = T_GS,
                 "p" = pchisq(T_GS, df = the_df, lower.tail = FALSE),
                 "mle0" = the_mle,
                 "mle1" = mle_alt,
                 "wald_cov" = robust_wald_cov,
                 "wald_se" = robust_wald_se,
                 "coef_tab" = coef_tab,
                 "penalty" = penalty)

  return(structure(result, class = "raoFit"))

}
