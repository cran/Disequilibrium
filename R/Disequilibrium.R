

#' TransformSigma_PDtoR3
#'
#' @param vec A length 3 vector of the variance of equation 1, followed by the
#' covariance of equations 1 and 2, followed by the variance of equation 2.
#'
#' @return A length 3 vector spanning unrestricted R3.
#'
#' @export
#'
#' @examples
#'PD_vec <- c(1, 0, 1)
#'TransformSigma_PDtoR3(PD_vec)
#'
TransformSigma_PDtoR3 <- function(vec){
  vec[2] <- atanh(vec[2] / (sqrt(vec[1]) * sqrt(vec[3])))
  vec[1] <- log(vec[1])
  vec[3] <- log(vec[3])
  return(vec)
}


#' TransformSigma_R3toPD
#'
#' @param vec A length 3 vector output from TransformSigma_PDtoR3.
#'
#' @return A length 3 vector spanning a positive definite space.
#'
#' @export
#'
#' @examples
#'PD_vec <- c(1, 0, 1)
#'R3_vec <- TransformSigma_PDtoR3(PD_vec)
#'TransformSigma_R3toPD(R3_vec)
#'
TransformSigma_R3toPD <- function(vec){
  vec[1] <- exp(vec[1])
  vec[3] <- exp(vec[3])
  vec[2] <- tanh(vec[2]) * sqrt(vec[1]) * sqrt(vec[3])
  return(vec)
}

#' Derivative of likelihood with respect to the inverse hyperbolic tangent of correlation
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar log of inverse hyperbolic tangent of the correlation of equations 1 and 2.
#'
#' @return A vector of derivatives for each observation.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'
#'d = DlhoodDatanhrho(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2])
#'head(d)
#'
DlhoodDatanhrho <- function(Y,mu,logsigma11,logsigma22,atanhrho){
  -stats::dnorm((Y-mu[,1] - tanh(atanhrho)*sqrt(exp(logsigma11))*(Y-mu[,2])/sqrt(exp(logsigma22)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma11))^2)) *
    stats::dnorm(Y,mu[,2],sqrt(exp(logsigma22))) *
    ((Y-mu[,1])*sinh(atanhrho)/sqrt(exp(logsigma11)) - (Y-mu[,2])*cosh(atanhrho)/sqrt(exp(logsigma22))) +
    -stats::dnorm((Y-mu[,2] - tanh(atanhrho)*sqrt(exp(logsigma22))*(Y-mu[,1])/sqrt(exp(logsigma11)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2)) *
    stats::dnorm(Y,mu[,1],sqrt(exp(logsigma11))) *
    ((Y-mu[,2])*sinh(atanhrho)/sqrt(exp(logsigma22)) - (Y-mu[,1])*cosh(atanhrho)/sqrt(exp(logsigma11)))
}


#' Derivative of log likelihood with respect to the inverse hyperbolic tangent of correlation
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar of the inverse hyperbolic tangent of the correlation of equations 1 and 2.
#' @param lhood A vector of length \eqn{N}{N} of likelihood values.
#'
#' @return A vector of derivatives for each observation.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'lhood = exp(-LLikelihoodDE(theta, Y, X, ensurePD = TRUE, summed = FALSE))
#'
#'d <- DllhoodDatanhrho(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2], lhood = lhood)
#'head(d)
#'
DllhoodDatanhrho <- function(Y,mu,logsigma11,logsigma22,atanhrho,lhood){

  (1/lhood) * DlhoodDatanhrho(Y,mu,logsigma11,logsigma22,atanhrho)

}

#' Derivative of likelihood with respect to the log of variance for equation 1
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar of the inverse hyperbolic tangent of the correlation of equations 1 and 2.
#'
#' @return A vector of derivatives for each observation.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'
#'d = DlhoodDlogsigma11(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2])
#'head(d)
#'
DlhoodDlogsigma11 <- function(Y,mu,logsigma11,logsigma22,atanhrho){
  .5*(-stats::dnorm((Y-mu[,1] - tanh(atanhrho)*sqrt(exp(logsigma11))*(Y-mu[,2])/sqrt(exp(logsigma22)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma11))^2)) *
        stats::dnorm(Y,mu[,2],sqrt(exp(logsigma22))) *
        -(Y-mu[,1])/sqrt(exp(logsigma11)*(1-tanh(atanhrho)^2)) +
        (-stats::dnorm((Y-mu[,2] - tanh(atanhrho)*sqrt(exp(logsigma22))*(Y-mu[,1])/sqrt(exp(logsigma11)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2)) *
           stats::dnorm(Y,mu[,1],sqrt(exp(logsigma11)))*
           tanh(atanhrho)*(Y-mu[,1])*sqrt(exp(logsigma22))/sqrt(exp(logsigma11))/sqrt(exp(logsigma22)*(1-tanh(atanhrho)^2))+
           (1-stats::pnorm(Y,mu[,2] + tanh(atanhrho)*sqrt(exp(logsigma22))*(Y-mu[,1])/sqrt(exp(logsigma11)),sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2)))*
           stats::dnorm(Y,mu[,1],sqrt(exp(logsigma11))) *
           (exp(-logsigma11)*(Y-mu[,1])^2-1)))
}

#' Derivative of log likelihood with respect to the log of variance for equation 1
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar of the inverse hyperbolic tangent of the correlation of equations 1 and 2.
#' @param lhood A vector of length \eqn{N}{N} of likelihood values.
#'
#' @return A vector of derivatives for each observation.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'lhood = exp(-LLikelihoodDE(theta, Y, X, ensurePD = TRUE, summed = FALSE))
#'
#'d <- DllhoodDlogsigma11(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2], lhood = lhood)
#'head(d)
#'
DllhoodDlogsigma11 <- function(Y,mu,logsigma11,logsigma22,atanhrho,lhood){

  (1/lhood) * DlhoodDlogsigma11(Y,mu,logsigma11,logsigma22,atanhrho)

}


#' Derivative of likelihood with respect to the coefficients of equation 1
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar of the inverse hyperbolic tangent of the correlation of equations 1 and 2.
#' @param X1 A \eqn{N \times k_1}{N x k[1]} design matrix for equation 1.
#'
#' @return A matrix of derivatives for each observation and parameter.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'
#'d = DlhoodDbeta1(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2], X1 = X1)
#'head(d)
#'
DlhoodDbeta1 <- function(Y,mu,logsigma11,logsigma22,atanhrho,X1){
  matrix(
    -stats::dnorm((Y-mu[,1] - tanh(atanhrho)*sqrt(exp(logsigma11))*(Y-mu[,2])/sqrt(exp(logsigma22)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma11))^2)) *
      stats::dnorm(Y,mu[,2],sqrt(exp(logsigma22))),nrow(X1),ncol(X1)) *
    -X1/sqrt(exp(logsigma11)*(1-tanh(atanhrho)^2)) +

    matrix(
      -stats::dnorm((Y-mu[,2] - tanh(atanhrho)*sqrt(exp(logsigma22))*(Y-mu[,1])/sqrt(exp(logsigma11)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2)) *
        stats::dnorm(Y,mu[,1],sqrt(exp(logsigma11)))*
        (tanh(atanhrho)*sqrt(exp(logsigma22))/sqrt(exp(logsigma11)))/sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2),nrow(X1),ncol(X1)) * X1 +

    matrix(
      (1-stats::pnorm(Y,mu[,2] + tanh(atanhrho)*sqrt(exp(logsigma22))*(Y-mu[,1])/sqrt(exp(logsigma11)),sqrt((1-tanh(atanhrho)^2)*sqrt(exp(logsigma22))^2)))*
        stats::dnorm(Y,mu[,1],sqrt(exp(logsigma11))) *
        exp(-logsigma11)*(Y-mu[,1]),nrow(X1),ncol(X1))*(X1)
}

#' Gradient of log likelihood with respect to the coefficients of equation 1
#'
#' @param Y A vector of observed responses.
#' @param mu A \eqn{N \times 2}{N x 2} matrix of means for equations 1 and 2.
#' @param logsigma11 A scalar log of the variance of the equation 1.
#' @param logsigma22 A scalar log of the variance of the equation 2.
#' @param atanhrho A scalar of the inverse hyperbolic tangent of the correlation of equations 1 and 2.
#' @param lhood A vector of length \eqn{N}{N} of likelihood values.
#' @param X1 A \eqn{N \times k_1}{N x k[1]} design matrix for equation 1.
#'
#' @return A matrix of derivatives for each observation and parameter.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
#'lhood = exp(-LLikelihoodDE(theta, Y, X, ensurePD = TRUE, summed = FALSE))
#'
#'d = DllhoodDbeta1(Y = Y, mu = mu, logsigma11 = theta[p1 + p2 + 1],
#'    logsigma22 = theta[p1 + p2 + 3], atanhrho = theta[p1 + p2 + 2], lhood = lhood, X1 = X1)
#'head(d)
#'
DllhoodDbeta1 <- function(Y,mu,logsigma11,logsigma22,atanhrho,lhood,X1){

  (1/matrix(lhood,nrow(X1),ncol(X1))) * DlhoodDbeta1(Y,mu,logsigma11,logsigma22,atanhrho,X1)

}

#' Gradient of log likelihood with respect to all parameters
#'
#' @param theta A vector of parameter values to obtain the gradient at. The
#' order of parameters is coefficients of equation 1, coefficients of equation 2,
#' log variance of equation 1, inverse hyperbolic tangent of the correlation of equations 1 and 2,
#' and log variance of equation 2.
#' @param Y A vector of observed responses.
#' @param X A list of two elements. The first element is a \eqn{N \times k_1}{N x k[1]}
#' design matrix for equation 1 and the second element is a \eqn{N \times k_2}{N x k[2]}
#' design matrix for equation 2.
#' @param summed A logical to determine if gradient values are summed over observations.
#'
#' @return A \eqn{(k_1 + k_2 + 3)}{(k[1] + k[2] + 3)} dimension vector of
#' derivatives if \code{summed = TRUE}, else a
#' \eqn{N \times (k_1 + k_2 + 3)}{N x (k[1] + k[2] + 3)} matrix of derivatives.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'
#'Gradient = GradientDE(theta, Y, X, summed = TRUE)
#'head(Gradient)
#'
GradientDE <- function(theta, Y, X, summed = TRUE){

  p1 = ncol(X[[1]])
  p2 = ncol(X[[2]])
  p = p1 + p2

  # Set up parameters
  mu = cbind(X[[1]] %*% theta[1:p1], X[[2]] %*% theta[(p1 + 1):(p1 + p2)])
  logsigma11 = theta[p+1]
  atanhrho = theta[p+2]
  logsigma22 = theta[p+3]

  # Precompute likelihood
  lhood = exp(-LLikelihoodDE(theta, Y, X, ensurePD = TRUE, summed = FALSE))

  gradvec = cbind(
    DllhoodDbeta1(Y,mu,logsigma11,logsigma22,atanhrho,lhood,X[[1]]),
    DllhoodDbeta1(Y,mu[,2:1],logsigma22,logsigma11,atanhrho,lhood,X[[2]]),
    DllhoodDlogsigma11(Y,mu,logsigma11,logsigma22,atanhrho,lhood),
    DllhoodDatanhrho(Y,mu,logsigma11,logsigma22,atanhrho,lhood),
    DllhoodDlogsigma11(Y,mu[,2:1],logsigma22,logsigma11,atanhrho,lhood))

  if(summed){
    return(colSums(gradvec))
  } else {
    return(gradvec)
  }

}


#' Negative gradient of log likelihood with respect to all parameters
#'
#' @param theta A vector of parameter values to obtain the gradient at. The
#' order of parameters is coefficients of equation 1, coefficients of equation 2,
#' log variance of equation 1, inverse hyperbolic tangent of the correlation of equations 1 and 2,
#' and log variance of equation 2.
#' @param Y A vector of observed responses.
#' @param X A list of two elements. The first element is a \eqn{N \times k_1}{N x k[1]}
#' design matrix for equation 1 and the second element is a \eqn{N \times k_2}{N x k[2]}
#' design matrix for equation 2.
#' @param summed A logical to determine if gradient values are summed over observations.
#'
#' @return A \eqn{(k_1 + k_2 + 3)}{(k[1] + k[2] + 3)} dimension vector of
#' derivatives if \code{summed = TRUE}, else a
#' \eqn{N \times (k_1 + k_2 + 3)}{N x (k[1] + k[2] + 3)} matrix of derivatives.
#'
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'
#'Gradient = nGradientDE(theta, Y, X, summed = TRUE)
#'head(Gradient)
#'
nGradientDE <- function(theta, Y, X, summed = TRUE){
  return(-GradientDE(theta, Y, X, summed = summed))
}


#' Negative log likelihood of market in disequilibrium model
#'
#' @param theta A vector of parameter values to obtain the gradient at. The
#' order of parameters is coefficients of equation 1, coefficients of equation 2,
#' variance of equation 1, correlation of equations 1 and 2,
#' and variance of equation 2.
#' @param Y A vector of observed responses.
#' @param X A list of two elements. The first element is a \eqn{N \times k_1}{N x k[1]}
#' design matrix for equation 1 and the second element is a \eqn{N \times k_2}{N x k[2]}
#' design matrix for equation 2.
#' @param ensurePD A logical  to determine if the covariance matrix is
#' transformed to an unrestricted 3 dimension real space (\code{ensurePD = TRUE})
#' or not (\code{ensurePD = FALSE}).
#' @param summed A logical to determine if the negative log likelihood values
#' are summed over observations.
#'
#' @return A scalar value of the negative log likelihood if \code{summed = TRUE},
#' else a \eqn{N}{N} length vector of negative log likelihood observations.
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'
#'p1 = 2
#'p2 = 2
#'theta = c(beta01, beta02, log(SigmaX[1, 1]), atanh(SigmaX[1, 2]), log(SigmaX[2, 2]))
#'
#'lhood = LLikelihoodDE(theta, Y, X, summed = TRUE)
#'head(LLikelihoodDE)
#'
LLikelihoodDE <- function(theta, Y, X, ensurePD = TRUE, summed = TRUE){
  p1 <- ncol(X[[1]])
  p2 <- ncol(X[[2]])
  p <- p1 + p2
  beta <- list(beta1 = theta[1:p1], beta2 = theta[(p1 + 1):p])
  Sigma <- matrix(c(theta[p + 1], theta[p + 2], theta[p + 2], theta[p + 3]), 2, 2)

  # convert covariance matrix to normal parameter space
  if(ensurePD){
    newSigma <- TransformSigma_R3toPD(c(Sigma[1, 1], Sigma[1, 2], Sigma[2, 2]))
    Sigma[1, 1] <- newSigma[1]
    Sigma[1, 2] <- newSigma[2]
    Sigma[2, 1] <- Sigma[1, 2]
    Sigma[2, 2] <- newSigma[3]
  }

  # precomputation
  XBeta1 <- X[[1]] %*% beta[[1]]
  XBeta2 <- X[[2]] %*% beta[[2]]

  # compute likelihood
  PrY <-
    (1 - stats::pnorm(q = Y,
                      mean = XBeta1 + Sigma[1, 2] * (Y - XBeta2) / Sigma[2, 2],
                      sd = sqrt(Sigma[1, 1] - (Sigma[1, 2] / sqrt(Sigma[2, 2])) ^ 2))) *
    stats::dnorm(x = Y,
                 mean = XBeta2,
                 sd = sqrt(Sigma[2, 2])) +
    (1 - stats::pnorm(q = Y,
                      mean = XBeta2 + Sigma[1, 2] * (Y - XBeta1) / Sigma[1, 1],
                      sd = sqrt(Sigma[2, 2] - (Sigma[1, 2] / sqrt(Sigma[1, 1])) ^ 2))) *
    stats::dnorm(Y, XBeta1, sqrt(Sigma[1, 1]))

  llhood <- log(PrY)
  if(summed){
    return(-sum(llhood))
  } else {
    return(-llhood)
  }
}


#' @title Market in Disequilibrium Model
#'
#' @description
#'
#' \code{DE} estimates a market in disequilibrium model.
#'
#' The market in disequilibrium model is defined as follows.
#' Let \eqn{i} denote the \eqn{i}th observation which takes values from \eqn{1}
#' to \eqn{N}, \eqn{X_1}{X[1]} be a covariate matrix of dimension
#' \eqn{N \times k_1}{N x k[1]}, \eqn{X_2}{X[2]} be a covariate matrix of
#' dimension \eqn{N \times k_2}{N x k[2]}, \eqn{X_{1i}}{X[1i]} be the \eqn{i}th
#' row of \eqn{X_1}{X[1]}, \eqn{X_{2i}}{X[2i]} be the \eqn{i}th row of
#' \eqn{X_2}{X[2]}, \eqn{\beta_1}{\beta[1]} be a coefficient vector of length
#' \eqn{k_1}{k[1]} and \eqn{\beta_2}{\beta[2]} be a coefficient vector of length
#' \eqn{k_2}{k[2]}. Define the latent response for stage one to be
#' \deqn{y_{1i}^\star = X_{1i} \beta_1 + \epsilon_{1i}}{y[1i]* = X[1i] \beta[1] + \epsilon[1i]}
#' and stage two to be
#' \deqn{y_{2i}^\star = X_{2i} \beta_2 + \epsilon_{2i}.}{y[2i]* = X[2i] \beta[2] + \epsilon[2i].}
#' Define the observed outcome to be
#' \eqn{y_{i}=min(y_{1i}^\star, y_{2i}^\star).}{y[i]=min(y[1i]*,y[2i]*).} The pair
#'  \eqn{(\epsilon_{1i},\epsilon_{2i})}{(\epsilon[1i],\epsilon[2i])} is distributed independently
#'  and identically multivariate normal with means
#'  \eqn{E[\epsilon_{1i}] = E[\epsilon_{2i}] = 0}{E(\epsilon[1i]) = E(\epsilon[2i]) = 0},
#'  variances
#'  \eqn{Var[\epsilon_{1i}] = \sigma_{11}}{Var(\epsilon[1i]) = \sigma[11]}, \eqn{Var[\epsilon_{2i}] = \sigma_{22}}{Var(\epsilon[2i]) = \sigma[22]},
#'  and covariance
#'  \eqn{Cov(\epsilon_{1i},\epsilon_{2i}) = \sigma_{12}}{Cov(\epsilon[1i],\epsilon[2i]) = \sigma[12]}.
#'
#'  The model is estimated by
#'  (frequentist) maximum likelihood. The default maximum likelihood algorithm is based off the Limited-memory
#'  Broyden Fletcher Goldfarb Shanno (L-BFGS) algorithm. See
#'  \link[optimr]{optimr} for details.
#'
#' @param formula An object of class \link[Formula]{Formula}: a symbolic
#'  description of the model to be fitted. The details of model specification
#'  are given under 'Details'.
#' @param data A required data frame containing the variables provided in \code{formula}.
#' @param subset An optional numeric vector specifying a subset of observations to be used in the fitting process.
#' @param par A vector of initial values for the parameters for which optimal values are to be found.
#'  The order of the parameters is the coefficients of equation 1, the coefficients of equation 2, the variance for equation 1,
#'  the covariance of equations 1 and 2, and the variance of equation 2.  The default is 0 for all parameters, except the variance
#'  parameters, which are set to 1.
#' @param control A list of control parameters. See 'Details'.
#'
#' @return
#'
#' An object of class 'DE' is returned as a list with the following components:
#'
#' \describe{
#'
#' \item{par}{The set of parameter maximum likelihood estimates.}
#'
#' \item{value}{The negative of the log likelihood evaluated at \code{par}.}
#'
#' \item{counts}{A two-element integer vector giving the number of calls to \code{fn} and \code{gr} respectively.
#'  This excludes those calls needed to compute the Hessian, if requested, and any calls to \code{fn} to compute a finite-difference
#'  approximation to the gradient.}
#'
#' \item{convergence}{An integer code. '0' indicates successful completion.}
#'
#' \item{message}{A character string giving any additional information returned by the optimizer, or \code{NULL}.}
#'
#' \item{Sigma}{The estimated covariance matrix.  This is a subset of \code{par}.}
#'
#' \item{Xbeta1}{The predicted outcome for each observation in equation 1 using the parameters estimated in \code{par}.}
#'
#' \item{Xbeta2}{The predicted outcome for each observation in equation 2 using the parameters estimated in \code{par}.}
#'
#' \item{hessian}{A numerical approximation to the Hessian matrix of the likelihood function at the estimated parameter values.
#'  The Hessian is returned even if it is not requested.}
#'
#' }
#'
#'
#' @export
#'
#' @details
#'
#' Models for \code{DE} are specified symbolically. A typical
#'  model has the form \code{response ~ terms1 | terms2} where \code{response}
#'  is the name of the numeric response variable and \code{terms1} and \code{terms2}
#'  are each a series of variables that specifies the linear predictor(s) of the respective model.
#'  For example, if the first equation has two independent variables, \code{X1} and \code{X2},
#'  then \code{terms1 = X1 + X2}.
#'
#' A \code{formula} has an implied intercept term for both equations. To remove the
#'  intercept from equation 1 use either \code{response ~ terms1 - 1 | terms2}
#'  or \code{response ~ 0 + terms1 | terms2}. The intercept may be removed from equation 2 analgously.
#'
#' The \code{control} argument is a list that can supply any of the following components:
#'
#' \describe{
#'
#' \item{method}{A method to be used in the function \code{optimr}. The default is \code{"L-BFGS-B"}.
#'  See \link[optimr]{optimr} for further details.}
#'
#' \item{lower, upper}{Bounds on the variables for methods such as \code{"L-BFGS-B"} that can handle box
#'  (or bounds) constraints. The default is \code{-Inf} and \code{Inf}, respectively.}
#'
#' \item{hessian}{A logical control that if TRUE forces the computation of an approximation to the Hessian
#'  at the final set of parameters. See \link[optimr]{optimr} for further details. The default is FALSE.}
#'
#' \item{na.action}{A function indicating what happens when data contains NAs. The default is \link{na.omit}.
#'  The only other possible value is \link{na.fail}.}
#'
#' \item{ensurePD}{A logical. If TRUE, the covariance matrix will be manually converted to a
#'  positive definite matrix for optimization. The default is TRUE.}
#'
#' \item{Equation1Name}{A string name for the first equation. The default is "_1".}
#'
#' \item{Equation2Name}{A string name for the second equation. The default is "_2".}
#'
#'}
#'
#' An object of class 'DE' is returned as a list with the following attributes:
#'
#' \describe{
#'
#' \item{originalNamesX1}{The original names of the covariates for equation 1 specified in \code{formula}.}
#'
#' \item{originalNamesX2}{The original names of the covariates for equation 2 specified in \code{formula}.}
#'
#' \item{namesX1}{The names of the covariates for equation 1 to be used in the \code{summary.DE} function.
#'  This is equivalent to \code{paste0(originalNamesX1, Equation1Name)}.}
#'
#' \item{namesX2}{The names of the covariates for equation 2 to be used in the \code{summary.DE} function.
#'  This is equivalent to \code{paste0(originalNamesX2, Equation2Name)}.}
#'
#' \item{namesSigma}{The names of the variance-covariance matrix parameters to be used in the \code{summary.DE} function.}
#'
#' \item{Equation1Name}{The user-specified name of equation 1.}
#'
#' \item{Equation2Name}{The user-specified name of equation 2.}
#'
#' \item{betaN1}{The number of coefficient and slope parameters to be estimated in equation 1.}
#'
#' \item{betaN2}{The number of coefficient and slope parameters to be estimated in equation 2.}
#'
#' }
#'
#'
#' @references
#' \cite{Gourieroux, C. (2000). Econometrics of Qualitative Dependent Variables (Themes in Modern Econometrics) (P. Klassen, Trans.). Cambridge: Cambridge University Press. doi:10.1017/CBO9780511805608}
#'
#' \cite{Maddala, G. (1983). Limited-Dependent and Qualitative Variables in Econometrics (Econometric Society Monographs). Cambridge: Cambridge University Press. doi:10.1017/CBO9780511810176}
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'df = data.frame(Y = Y, X1 = X1, X2 = X2)
#'
#'results = DE(formula = Y ~ X1 | X2, data = df)
#'
#'str(results)
#'
DE <- function(formula, data, subset = NULL, par = NULL, control = list()){
  # prep inputs
  mf <- match.call(expand.dots = FALSE)
  mf <- mf[c(1, match(c('formula', 'data', 'subset', 'na.action'), names(mf), 0))]
  f <- Formula::Formula(formula)
  mf[[1]] <- as.name('model.frame')
  mf$formula <- f
  mf <- eval(mf, parent.frame())
  if(!is.null(subset)){
    data <- data[subset,]
  }
  if(is.null(control$na.action)){
    control$na.action <- 'na.omit'
  }
  if(control$na.action == 'na.omit'){
    data <- stats::na.omit(data)
  }
  if(control$na.action == 'na.fail'){
    data <- stats::na.fail(data)
  }
  if(is.null(control$method)){
    control$method <- 'L-BFGS-B'
  }
  if(is.null(control$lower)){
    control$lower <- -Inf
  }
  if(is.null(control$upper)){
    control$upper <- Inf
  }
  if(is.null(control$hessian)){
    control$hessian <- F
  }
  Y <- stats::model.response(mf)
  X1 <- stats::model.matrix(f, data = mf, rhs = 1)
  X2 <- stats::model.matrix(f, data = mf, rhs = 2)
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  p <- p1 + p2
  originalNamesX1 <- colnames(X1)
  originalNamesX2 <- colnames(X2)

  # check inputs
  if(length(Y) == 0 | nrow(X1) == 0 | nrow(X2) == 0){
    stop('formula must be of the form Y ~ X1 | X2')
  }
  if(!is.numeric(c(Y, X1, X2))){
    stop('variables in formula must be numeric')
  }
  if(!is.null(par)){
    if(length(par) != p + 3){
      stop(paste0('par must be a vector of length equal to the number of free parameters (', p + 3, ')'))
    }
  }
  if(!is.null(par) & !is.numeric(par)){
    stop('par must be a numeric vector')
  }
  if(!is.null(par)){
    if((par[p + 1] <= 0) |
       (par[p + 3] <= 0) |
       (par[p + 3] - par[p + 2] / par[p + 1] < 0)){
      stop('variance-covariance parameters in par must correspond to a positive definite matrix')
    }
  }
  if(!is.numeric(control$lower)){
    stop('control$lower must be numeric')
  }
  if(!is.numeric(control$upper)){
    stop('control$upper must be numeric')
  }
  if(!is.logical(control$hessian)){
    stop('control$hessian must be a logical')
  }
  if(!(control$na.action %in% c('na.omit', 'na.fail'))){
    stop('control$na.action must be "na.omit" or "na.fail"')
  }
  if(is.null(control$Equation1Name)){
    control$Equation1Name <- '_1'
  }
  if(is.null(control$Equation2Name)){
    control$Equation2Name <- '_2'
  }

  # prep data
  if(is.null(par)){
    par <- c(rep(0, p), 1, 0, 1)
  }
  par[(p + 1):(p + 3)] <- TransformSigma_PDtoR3(par[(p + 1):(p + 3)])

  # MLE
  out <- optimr::optimr(par = par, fn = LLikelihoodDE, gr = nGradientDE, lower = control$lower,
                        upper = control$upper,
                        method = control$method, hessian = control$hessian,
                        X = list(X1, X2), Y = Y, control = control)
  out$par[(p + 1):(p + 3)] <- TransformSigma_R3toPD(out$par[(p + 1):(p + 3)])

  # create DE class
  class(out) <- c(class(out), 'DE')
  attr(out, 'originalNamesX1') <- originalNamesX1
  attr(out, 'originalNamesX2') <- originalNamesX2
  attr(out, 'namesX1') <- paste0(originalNamesX1, control$Equation1Name)
  attr(out, 'namesX2') <- paste0(originalNamesX2, control$Equation2Name)
  attr(out, 'namesSigma') <- c(paste0('Variance', control$Equation1Name), 'Covariance', paste0('Variance', control$Equation2Name))
  attr(out, 'Equation1Name') <- control$Equation1Name
  attr(out, 'Equation2Name') <- control$Equation2Name
  attr(out, 'betaN1') <- p1
  attr(out, 'betaN2') <- p2
  out$Sigma <- matrix(c(out$par[p + 1], out$par[p + 2], out$par[p + 2], out$par[p + 3]), nrow = 2, ncol = 2)
  Xbeta1 <- X1 %*% matrix(out$par[1:p1], nrow = p1, ncol = 1)
  attributes(Xbeta1)$dimnames <- NULL
  out$Xbeta1 <- Xbeta1
  Xbeta2 <- X2 %*% matrix(out$par[(p1 + 1):(p1 + p2)], nrow = p2, ncol = 1)
  attributes(Xbeta2)$dimnames <- NULL
  out$Xbeta2 <- Xbeta2
  out$hessian <- numDeriv::hessian(func = LLikelihoodDE, x = out$par,
                                   X = list(X1, X2), Y = Y, ensurePD = FALSE)
  return(out)
}

#' @title Summary method for class 'DE'
#'
#' @param object An object of class \code{DE} for which a summary is desired.
#' @param ... Unused
#'
#' @return
#'
#' A matrix summary of estimates is returned.  The columns are:
#'
#' \describe{
#'
#' \item{Estimate}{Maximum likelihood point estimate.}
#'
#' \item{Std. Error}{Asymptotic standard error estimate of maximum likelihood point estimators using numerical hessian.}
#'
#' \item{z value}{z value for zero value null hypothesis using asymptotic standard error estimate.}
#'
#' \item{Pr(>|z|)}{P value for a two sided null hyptothesis test using the z value.}
#'
#' }
#'
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'df = data.frame(Y = Y, X1 = X1, X2 = X2)
#'
#'results = DE(formula = Y ~ X1 | X2, data = df)
#'
#'summary(results)
#'
summary.DE <- function(object,...){
  dots <- match.call(expand.dots = TRUE)
  if(is.null(dots$digits)){
    dots$digits <- max(3, getOption("digits") - 3)
  }
  coeff <- object$par
  if(!is.null(object$hessian)){
    se <- sqrt(diag(solve(object$hessian)))
  } else {
    se <- rep(NA,length(coeff))
  }
  zstat <- abs(coeff / se)
  pval <- stats::pnorm(zstat, lower.tail = FALSE) * 2
  out <- cbind(coeff, se, zstat, pval)
  colnames(out) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  rownames(out) <- c(attr(object, 'namesX1'), attr(object, 'namesX2'), attr(object, 'namesSigma'))
  return(out)
}

#' @title Predict method for class 'DE'
#'
#' @param object An object of class \code{DE}.
#' @param newdata An optional data frame with column names matching the dependent variables specified in
#'  the \code{formula} of the \code{DE} function. If not provided, the \code{data} from the \code{DE} function
#'  will be used.
#' @param ... Unused
#'
#' @return
#' A data frame is returned.  The columns are:
#'
#' \describe{
#'
#' \item{Y_1}{Linear prediction of the outcome variable in equation 1.}
#'
#' \item{Y_2}{Linear prediction of the outcome variable in equation 2.}
#'
#' \item{Min(Y_1,Y_2)}{The minimum of \code{Y_1} and \code{Y_2}.}
#'
#' \item{Prob(Y_2>Y_1)}{The probability that the outcome variable in equation 2 is greater than the outcome
#' variable in equation 1.  This is the probability that \code{Y_1} is the observed quantity.}
#'
#' }
#'
#' @export
#'
#' @examples
#'set.seed(1775)
#'library(MASS)
#'beta01 = c(1,1)
#'beta02 = c(-1,-1)
#'N = 10000
#'SigmaEps = diag(2)
#'SigmaX = diag(2)
#'MuX = c(0,0)
#'par0 = c(beta01, beta02, SigmaX[1, 1], SigmaX[1, 2], SigmaX[2, 2])
#'
#'Xgen = mvrnorm(N,MuX,SigmaX)
#'X1 = cbind(1,Xgen[,1])
#'X2 = cbind(1,Xgen[,2])
#'X = list(X1 = X1,X2 = X2)
#'eps = mvrnorm(N,c(0,0),SigmaEps)
#'eps1 = eps[,1]
#'eps2 = eps[,2]
#'Y1 = X1 %*% beta01 + eps1
#'Y2 = X2 %*% beta02 + eps2
#'Y = pmin(Y1,Y2)
#'df = data.frame(Y = Y, X1 = X1, X2 = X2)
#'
#'results = DE(formula = Y ~ X1 | X2, data = df)
#'
#'head(predict(results))
#'
predict.DE <- function(object, newdata = NULL,...){
  if(is.null(newdata)){
    Xbeta1 <- object$Xbeta1
    Xbeta2 <- object$Xbeta2
  } else {
    X1 <- matrix(newdata[, attr(object, 'originalNamesX1')])
    X2 <- matrix(newdata[, attr(object, 'originalNamesX2')])
    Xbeta1 <- X1 %*% matrix(object$par[1:attr(object, 'betaN1')], nrow = attr(object, 'betaN1'), ncol = 1)
    Xbeta2 <- X2 %*% matrix(object$par[(attr(object, 'betaN1') + 1):(attr(object, 'betaN1') + attr(object, 'betaN2'))],
                            nrow = attr(object, 'betaN1'), ncol = 1)
  }
  out <- data.frame(predictedY1 = Xbeta1,
                    predictedY2 = Xbeta2,
                    predictedMinY = pmin(Xbeta1, Xbeta2),
                    prob_Y2_gt_Y1 = 1 - stats::pnorm(0,
                                                     mean = as.vector(c(-1, 1) %*% t(matrix(c(Xbeta1, Xbeta2), ncol = 2))),
                                                     sd = rep(sqrt(matrix(c(-1, 1), nrow = 1) %*%
                                                                     object$Sigma %*%
                                                                     matrix(c(-1, 1), ncol = 1)),
                                                              length(Xbeta1))))
  colnames(out) <- c(paste0('Y', attr(object, 'Equation1Name')),
                     paste0('Y', attr(object, 'Equation2Name')),
                     paste0('Min(Y', attr(object, 'Equation1Name'), ',Y', attr(object, 'Equation2Name'), ')'),
                     paste0('Prob(Y', attr(object, 'Equation2Name'), '>Y', attr(object, 'Equation1Name'), ')'))
  return(out)
}


