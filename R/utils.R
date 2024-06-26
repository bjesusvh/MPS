## Helping functions

## Calculate posterior BV in the form of mcmc samples.
## X A matrix of predictors of dimension n x k .
## B An array of mcmc samples for regression coeficients, M x k (predictors) x t.
## B0 A matrix of mcmc samples of M x t for intercept term in linear model.
## return An array of BVs in form of mcmc samples of dimension n x t x M
posteriorBV <- function(Xcand, B, B0){
  n <- nrow(Xcand)
  M <- dim(B)[1]
  t <- dim(B)[3]
  XB <- array(numeric(),c(n, t, M))
  for(m in 1:M){
    XB[ , , m ] <- Xcand %*% B[m , , ]
  }
  postBV <- sweep(XB, c(2,3), t(B0), "+")
  return(postBV)
}

# Function to simulate from multivariate normal distribution
# https://gallery.rcpp.org/articles/simulate-multivariate-normal/
mvrnormR <- function(n, mu, sigma) {
    ncols <- ncol(sigma)
    mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
    mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}

## Transform to a positive-definite matrix
toNearPD <- function(iter, x){
  mo <- xpnd(x[iter,])
  ms <- (mo + t(mo))/2
  m <- as.matrix(nearPD(ms)$mat)
  return(m)
}

## Loss functions
# Multivariate loss functions
multi_KL <- function(lower, upper, mu1, mu2, muS, K, Kinv){

  Z <- as.numeric(pmvnorm(lower = lower, upper = upper, mean = mu1, sigma = K)) # mvtnorm::

  if(Z < 1e-323){
    logZ = log(1e-323)
  }else{
      logZ = log(Z)
  }

  yppdf <- as.matrix(t(apply(mu2, 1, function(ypred) {mvrnormR(1, mu = ypred, sigma = K)})))

  S <- muS - mu1
  SKS <- as.numeric( t(S) %*% Kinv %*% S)   # S Pinv S
  muSmu2 <-  sweep(yppdf, 2, muS, '-')      # muS - mu2
  loss <- as.vector(apply(muSmu2, 1, function(x) {0.5*(t(x) %*% Kinv %*% x - SKS)}) - logZ)

  return(loss)
}

EnergyScore <- function(mu2, muS, K){
  yppdf <- as.matrix(t(apply(mu2, 1, function(ypred) {mvrnormR(1, mu = ypred, sigma = K)})))
  yppdf_p <- as.matrix(t(apply(mu2, 1, function(ypred) {mvrnormR(1, mu = ypred, sigma = K)})))
  muSmu2 <- sweep(yppdf, 2, muS, '-')     # mu2 - muS    or X-y
  mu2mu2 <- yppdf-yppdf_p                 # mu2 - mu2'   or X-X'
  xj_y <- as.vector(apply(muSmu2, 1, function(x)  sqrt(sum(x^2)) ))
  xj_xjp <- as.vector(apply(mu2mu2,1,function(x)  0.5*sqrt(sum(x^2))))
  loss <- xj_y - xj_xjp
  return(loss)
}

MALF <- function(mu2, muS, K, tau){
  yppdf <- as.matrix(t(apply(mu2, 1, function(ypred) {mvrnormR(1, mu = ypred, sigma = K)})))
  error <- -1*sweep(yppdf, 2, muS, '-')   # mu2 - muS
  abs.error <- abs(error)                 # abs del error
  sum.abs.e <- rowSums(abs.error)         # sum errors by rows
  loss <- sum.abs.e + as.vector(error %*% tau)
  return(loss)
}

# Aproximation of multivariate loss functions.
aproxMultiKL <- function(lower, upper, mu1, mu2, muS, K, Kinv){
  Z <- as.numeric(pmvnorm(lower = lower, upper = upper, mean = mu1, sigma = K))
  if(Z < 1e-323){
    logZ = log(1e-323)
  }else{
    logZ = log(Z)
  }
  S <- muS - mu1
  SKS <- as.numeric( t(S) %*% Kinv %*% S)
  muSmu2 <-  sweep(mu2, 2, muS, '-')  # muS - mu2
  loss <- as.vector(apply(muSmu2, 1, function(x) {0.5*(t(x) %*% Kinv %*% x - SKS)}) - logZ)
}

aproxEnergyScore <- function(mu2, muS){
  muSmu2 <- sweep(mu2, 2, muS, '-')       # mu2 - muS    or X-y
  loss <- as.vector(apply(muSmu2, 1, function(x)  sqrt(sum(x^2)) ))
}

aproxMALF <- function(mu2, muS, tau){
  error <- -1*sweep(mu2, 2, muS, '-')   # mu2 - muS
  abs.error <- abs(error)                 # abs del error
  sum.abs.e <- rowSums(abs.error)         # sum errors by rows
  loss <- sum.abs.e + as.vector(error %*% tau)
}

# https://stackoverflow.com/questions/37613345/r-convert-upper-triangular-part-of-a-matrix-to-symmetric-matrix
constuyeSimOnes <- function(m){
  m = m + t(m) - 2*diag(diag(m)) + diag(1,nrow=dim(m)[1])
  return (m)}

# # Univariate loss functions
# single_KL <- function(){
#
# }
#
# crps <- function(){
#
# }
#
# linlin <- function(){
#
# }





