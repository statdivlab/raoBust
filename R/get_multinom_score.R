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

    #set up betas
    betanonj <- rep(5, (p+1)*(J-2)+1) #vector of non \beta_j parameters, include all \beta_{k0} values and \beta_{k} for k \neq j, and \beta_{j0} as final element
    betaj <- rep(0, p) #\beta_j null value

    optim(par = betanonj, fn = multinom_mle_weak_null,, Y = Y, X = X, j = j)

    betanonj_null1mle <- tryCatch({optim(par = betanonj, fn = multinom_mle_weak_null,, Y = Y, X = X, j = j)$par},
                                    error = function(cond) {return(NA)}) #get optimal values of beta

    if (any(is.na(betanonj_null1mle))) {
      stop("Could not find MLE under H0")
    }

    beta_null1mle <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
    for (k in 1:(J-1)) {
      if (k < j) {
        beta_null1mle[ ,k] <- betanonj_null1mle[((k-1)*(p+1)+1):(k*(p+1))] #populate column k where k < j
      } else if (k == j) {
        beta_null1mle[,k] <- c(betanonj_null1mle[(p+1)*(J-2)+1], betaj) #populate j-th column
      } else if (k > j) {
        beta_null1mle[ ,k] <- betanonj_null1mle[((k-2)*(p+1)+1):((k-1)*(p+1))] #populate column k where k > j (have to shift further down by 1 in index since j not in vector)
      }
    }

    #jacobian of function of parameter h(\beta) = 0
    H1 <- matrix(data = rep(0, p*(p+1)*(J-1)), ncol = (p+1)*(J-1), nrow = p)
    H1[1:p, ((j-1)*(p+1) + 2):(j*(p+1))] <- diag(nrow = p, ncol = p)

    #terms necessary for computation of S, I, D matrices
    XaugBeta1 <- Xaug %*% beta_null1mle
    pJ <- (1 + rowSums(exp(XaugBeta1)))^(-1)
    ps <- as.vector(pJ)*exp(XaugBeta1)
    ps_full <- cbind(ps, pJ)


    #score of \beta evaluated at mle under null constraint
    S1 <- matrix(data = rep(0, (p+1)*(J-1)), ncol = 1) #initialize matrix for score
    for (k in 1:(J-1)) {
      S1[((k-1)*(p+1) + 1):(k*(p+1)) ,] <- colSums(as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug)
    }

    #information matrix evaluated at mle under null constraint
    I1 <- matrix(data = rep(0, (p+1)^2*(J-1)^2), ncol = (p+1)*(J-1)) #initialize matrix for information matrix
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        if (l < k || l > k) {
          I1[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <-  (as.vector(-ps_full[ ,k]*ps_full[ ,l]*N)*t(Xaug)) %*% Xaug
        } else if (l == k) {
          I1[((k-1)*(p+1) + 1):(k*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- (as.vector(-ps_full[ ,k]*(-1 + ps_full[ ,l])*N)*t(Xaug)) %*% Xaug
        }
      }
    }

    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D1 <- matrix(data = rep(0, length(S1)^2), ncol = length(S1))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S1_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S1_l <- as.vector(Y[ ,l] - ps_full[ ,l]*N)*Xaug
        D1[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S1_k)%*%S1_l
      }
    }

    the_mle <- betanonj_null1mle
    the_df <- p
    my_misc <- beta_null1mle

    #compute statistic!
    T_GS<- tryCatch({as.numeric(t(S1) %*% solve(I1) %*% t(H1) %*%
                                  (solve(H1 %*% solve(I1) %*% D1 %*% solve(I1) %*% t(H1))) %*%
                                  H1 %*% solve(I1) %*% S1)},
                             error = function(cond) {return(NA)})

  } else {
    #initialize value for beta_k0 for all k = 1, \dots, J-1
    betanots <- rep(1, J-1)

    #jacobian of function of parameter h(\beta) = 0
    H2 <- matrix(data = rep(0, p*(J-1)*(p+1)*(J-1)), ncol = (p+1)*(J-1), nrow = p*(J-1))
    for (k in 1:(J-1)) {
      H2[ ((k-1)*p + 1):(k*p), ((k-1)*(p+1) + 2):(k*(p+1))] <- diag(nrow = p, ncol = p)
    }

    #compute mle under null constraint
    betanots_null2mle <- optim(betanots, multinom_mle_strong_null, Y = Y, X = X)$par #get optimal values of beta_k0's

    ## AW TODO below not needed?
    # beta_null2mle <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
    # beta_null2mle[1, ] <- betanots_null2mle


    #terms necessary for computation of S, I, D matrices
    pJ <- (1 + sum(exp(betanots_null2mle)))^(-1)
    ps <- pJ*exp(betanots_null2mle)
    ps_full <- matrix(rep(c(ps, pJ), n), nrow = n, byrow = TRUE)


    #score of \beta evaluated at mle under null constraint
    S2 <- matrix(data = rep(0, (p+1)*(J-1)), ncol = 1) #initialize matrix for score
    for (k in 1:(J-1)) {
      S2[((k-1)*(p+1) + 1):(k*(p+1)) ,] <- colSums(as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug)
    }

    #information matrix evaluated at mle under null constraint
    I2 <- matrix(data = rep(0, (p+1)^2*(J-1)^2), ncol = (p+1)*(J-1)) #initialize matrix for information matrix
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        if (l < k || l > k) {
          I2[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (-as.vector(ps_full[ ,k]*ps_full[ ,l]*N)*Xaug)
        } else if (l == k) {
          I2[((k-1)*(p+1) + 1):(k*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (-as.vector(ps_full[ ,k]*(-1 + ps_full[ ,l])*N)*Xaug)
        }
      }
    }

    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D2 <- matrix(data = rep(0, length(S2)^2), ncol = length(S2))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S2_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S2_l <- as.vector(Y[ ,l] - ps_full[ ,l]*N)*Xaug
        D2[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S2_k)%*%S2_l
      }
    }

    the_mle <- betanots_null2mle
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
