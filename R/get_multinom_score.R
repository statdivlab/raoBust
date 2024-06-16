#' Robust score (Rao) tests for multinomial regression.
#'
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param strong This is by default specified as FALSE to compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0}.
#' If specified to be TRUE, the function instead computes the robust score statistic to test the strong null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0}.
#' @param j If `strong` is specified as FALSE, this argument must be supplied. This specifies for which category \eqn{j} you want to test the weak null hypothesis that \eqn{\beta_j = 0}.
#'
#' @return The robust score test statistic for the specified hypothesis test according to the strong and j parameters. TODO
#'
#' @importFrom stats nlm optim
#'
#' @author Shirley Mathur
#'
#' @export
get_multinom_score <- function(X, Y, strong = FALSE, j = NULL) {

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
      stop("Marginal test specified by user, but no argument to j provided to test null hypothesis of \beta_j = 0 for a user-specified category j.")
    }

    #get beta mle under weak null constraint, first setting up initial value of beta and then optimizing
    beta_init <- matrix(1, nrow = p + 1, ncol = J-1)
    beta_init[(2:(p+1)), j] <- 0
    beta_null1mle <- multinom_fisher_scoring(beta_init, X, Y, strong = FALSE, null_j = j)

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
    my_misc <- beta_null1mle

    #compute statistic!
    T_GS<- tryCatch({as.numeric(t(S1) %*% solve(I1) %*% t(H1) %*%
                                  (solve(H1 %*% solve(I1) %*% D1 %*% solve(I1) %*% t(H1))) %*%
                                  H1 %*% solve(I1) %*% S1)},
                             error = function(cond) {return(NA)})

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
    beta_null2mle <- multinom_fisher_scoring(beta_init, X, Y, strong = TRUE)

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
    my_misc <- NA # beta_null2mle

    #compute statistic!
    T_GS <- tryCatch({as.numeric(t(S2) %*% solve(I2) %*% t(H2) %*%
                                   (solve(H2 %*% solve(I2) %*% D2 %*% solve(I2) %*% t(H2))) %*%
                                   H2 %*% solve(I2) %*% S2)},
                     error = function(cond) {return(NA)})


  }
  
  
  #compute mle under alternative
  beta_alt <- rep(-0.02, (p+1)*(J-1))
  mle_alt <- matrix(optim(beta_alt, multinom_mle_alternative, Y = Y, X = X)$par,  nrow = p + 1, ncol = J-1)

  return(list("test_stat" = T_GS,
              "p" = pchisq(T_GS, df = the_df, lower.tail=FALSE),
              "mle0" = the_mle,
              "mle1" = mle_alt,
              "misc" = my_misc))

}
