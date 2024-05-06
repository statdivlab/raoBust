

#' Objective function to find MLE under weak null (\eqn{\beta_j = 0} for specific j)
#'
#' @param betanonj A vector containing the initial values for all \eqn{\beta_k}, with \eqn{k \neq j}, as well as all \eqn{\beta_{k0}, for k = 1, \dots, J-1}. 
#' In particular, this vector should be so that the first \eqn{(J-2)(p+1)} entries are \eqn{\beta_{10}, \beta_{1}^{\top}, \beta_{20}, \beta_{2}^{\top}, \dots, \beta_{(j-1)0}, \beta_{j-1}^{\top},  \beta_{(j+1)0}, \beta_{j+1}^{\top}, \dots,  \beta_{(J-1)0}, \beta_{J-1}^{\top}, \beta_{j0}}.
#' @param betaj This is the null hypothesized value for \eqn{\beta_j}, which is by default set to be \eqn{\beta_j = 0}.
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param j This specifies for which category you want to compute the MLE under the constraint that \eqn{\beta_j = 0}.
#' @return A vector containing the optimal betanonj values to maximize the log-likelihood under the null constraint that \eqn{\beta_j = 0}. The components are listed out in the same manner as in the betanonj parameter.
null1objective <- function(betanonj, betaj = rep(0, p), Y, X, j) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  beta <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
  for (k in 1:(J-1)) {
    if (k < j) {
      beta[ ,k] <- betanonj[((k-1)*(p+1)+1):(k*(p+1))]
    } else if (k == j) {
      beta[,k] <- c(betanonj[(p+1)*(J-2)+1], betaj)
    } else if (k > j) {
      beta[ ,k] <- betanonj[((k-2)*(p+1)+1):((k-1)*(p+1))]
    }
  }
  
  Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
  XaugBeta <- Xaug %*% beta
  pJ <- (1 + rowSums(exp(XaugBeta)))^(-1)
  ps <- as.vector(pJ)*exp(XaugBeta)
  ps_full <- cbind(ps, pJ)
  
  loglik <- sum(Y*log(ps_full))
  
  return(-loglik)
}


#' Objective function to find MLE under global null (\eqn{\beta_1 = beta_2 = \dots = \beta_{J-1} = 0)}
#'
#' @param betanots A vector containing the initial values for all \eqn{\beta_{k0}, for k = 1, \dots, J-1}. 
#' In particular, this vector should be so that the entries are ordered as \eqn{\beta_{10}, \beta_{20}, \dots, \beta_{(J-1)0}}.
#' @param Y This should be the n x J data matrix of outcomes.
#' @param X This should be the n x p design matrix of covariates.
#' @return A vector containing the optimal values for betanots  to maximize the log-likelihood under the null constraint that \eqn{\beta_1 = beta_2 = \dots = \beta_{J-1} = 0}. 
#' The components are listed out in the same manner as in the betanots parameter.

null2objective <- function(betanots, Y, X) {
  n <- nrow(Y) #get sample size
  J <- ncol(Y) #get number of taxa
  N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
  
  pJ <- (1 + sum(exp(betanots)))^(-1)
  ps <- pJ*exp(betanots)
  ps_full <- matrix(rep(c(ps, pJ), n), nrow = n, byrow = TRUE)
  
  loglik <- sum(Y*log(ps_full))
  
  
  return(-loglik)
}

#' Robust score (Rao) test for multinomial regression.
#'
#' @param Y This should be the \eqn{n x J} data matrix of outcomes.
#' @param X This should be the \eqn{n x p} design matrix of covariates.
#' @param joint This is by default specified as FALSE to compute the robust score statistic to test the weak null that for one specific \eqn{j}, \eqn{\beta_j = 0}. 
#' If specified to be TRUE, the function instead computes the robust score statistic to test the global null that \eqn{\beta_1 = \beta_2 = \dots = \beta_{J-1} = 0}.
#' @param j If `join` is specified as FALSE, this argument must be supplied. This specifies for which category \eqn{j} you want to test the weak null hypothesis that \eqn{\beta_j = 0}.
#' @return The robust score test statistic for the specified hypothesis test according to thejoint and j parameters.
get_multinom_score <- function(X, Y, joint = FALSE, j = NULL) {
  
  #get n, p, J values (used throughout rest of the function to compute relevant quantities)
  n <- nrow(Y)
  p <- ncol(X)
  J <- ncol(Y)
  
  #compute test statistics under weak null that \beta_j = 0 for user-specified j
  if (joint == FALSE) {
    #report error if user specified marginal test but does not supply an argument for j
    if(is.null(j)) {
      stop("Marginal test specified by user, but no argument to j provided to test null hypothesis of \beta_j = 0 for a user-specified category j.")
    }
    
    #set up betas
    betanonj <- rep(0, (p+1)*(J-2)+1) #vector of non \beta_j parameters, include all \beta_{k0} values and \beta_{k} for k \neq j, and \beta_{j0} as final element
    betaj <- rep(0, p) #\beta_j null value
      
    #jacobian of function of parameter h(\beta) = 0
    H1 <- matrix(data = rep(0, p*(p+1)*(J-1)), ncol = (p+1)*(J-1), nrow = p)
    H1[1:p, ((j-1)*(p+1) + 2):(j*(p+1))] <- diag(nrow = p, ncol = p)
      
      
    betanonj_null1mle <- tryCatch({nlm(f = null1objective, p = betanonj,Y = Y, X = X, j = j)$estimate},
                                    error = function(cond) {return(NA)}) #get optimal values of beta
      
    #return NA when optimization does not converge under null constraint
    if (any(is.na(betanonj_null1mle))) {
      return(NA)
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
      
    #terms necessary for computation of S, I, D matrices
    n <- nrow(Y)
    N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
    Xaug <- cbind(matrix(rep(1, n), ncol = 1), X) #augmented covariate matrix
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
          I1[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (as.vector(ps_full[ ,k]*ps_full[ ,l]*N)*Xaug)
        } else if (l == k) {
          I1[((k-1)*(p+1) + 1):(k*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (as.vector(ps_full[ ,k]*(-1 + ps_full[ ,l])*N)*Xaug)
        }
      }
    }
      
    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D1 <- matrix(data = rep(0, length(S1)^2), ncol = length(S1))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S1_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S1_l <- as.vector(Y[ ,l] - ps_full[ ,k]*N)*Xaug
        D1[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S1_k)%*%S1_l
      }
    }
      
    #compute statistic!
    T_GS<- tryCatch({as.numeric(t(S1) %*% solve(I1) %*% t(H1) %*% (solve(H1 %*% solve(I1) %*% D1 %*% solve(I1) %*% t(H1))) %*% H1 %*% solve(I1) %*% S1)},
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
    betanots_null2mle <- optim(betanots, null2objective, Y = Y, X = X)$par #get optimal values of beta_k0's
    beta_null2mle <- matrix(rep(0, (p+1)*(J-1)), ncol = J-1) #initialize full (p+1) x (J-1) beta matrix
    beta_null2mle[1, ] <- betanots_null2mle
    
    
    #terms necessary for computation of S, I, D matrices
    n <- nrow(Y)
    N <- matrix(apply(Y, MARGIN = 1, FUN = sum), nrow = n) #totals by sample
    Xaug <- cbind(rep(1,n), X)
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
          I2[((l-1)*(p+1) + 1):(l*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (as.vector(ps_full[ ,k]*ps_full[ ,l]*N)*Xaug)
        } else if (l == k) {
          I2[((k-1)*(p+1) + 1):(k*(p+1)),((k-1)*(p+1) + 1):(k*(p+1))] <- t(Xaug) %*% (as.vector(ps_full[ ,k]*(-1 + ps_full[ ,l])*N)*Xaug)
        }
      }
    }
    
    #D matrix (sum of S_iS_i^T for all i = 1, \dots, n)
    D2 <- matrix(data = rep(0, length(S2)^2), ncol = length(S2))
    for (k in 1:(J-1)) {
      for (l in 1:(J-1)) {
        S2_k <- as.vector(Y[ ,k] - ps_full[ ,k]*N)*Xaug
        S2_l <- as.vector(Y[ ,l] - ps_full[ ,k]*N)*Xaug
        D2[((k-1)*(p+1) + 1):(k*(p+1)) ,((l-1)*(p+1) + 1):(l*(p+1))] <- t(S2_k)%*%S2_l
      }
    }
    
    #compute statistic!
    T_GS <- tryCatch({as.numeric(t(S2) %*% solve(I2) %*% t(H2) %*% (solve(H2 %*% solve(I2) %*% D2 %*% solve(I2) %*% t(H2))) %*% H2 %*% solve(I2) %*% S2)},error = function(cond) {return(NA)})
    
    
  }
  
  return(T_GS)
  
}
