##############################################################################
##############################################################################
##############################################################################

## Quantitative Trait Selection

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



##############################################################################
##############################################################################
##############################################################################

## Ordinal Trait Selection

## Helping functions


# Función para calcular la divergencia KL entre dos distribuciones
compute_kl <- function(p, q) {
  if (length(p) != length(q)) stop("Las distribuciones p y q deben tener la misma longitud.")
  if (any(q == 0)) stop("La distribución de referencia q contiene ceros, lo cual no es válido para KL.")
  
  kl <- sum(p * log(p / q), na.rm = TRUE) # Calcula la divergencia KL
  return(kl)
}


# Función para calcular la KL promedio para un solo rasgo ordinal
compute_kl_single_trait <- function(mcmc_array, q) {
  # Verifica que la referencia q esté normalizada
  q <- q / sum(q) # Normaliza q si no lo está
  
  # Dimensiones del arreglo
  n <- dim(mcmc_array)[1] # Número de líneas
  k <- dim(mcmc_array)[2] # Número de categorías
  M <- dim(mcmc_array)[3] # Número de muestras MCMC
  
  if (length(q) != k) stop("La longitud de q debe coincidir con el número de categorías k.")
  
  # Calcula la KL promedio para cada línea
  kl_values <- numeric(n) # Vector para almacenar los valores KL
  for (i in 1:n) {
    kl_samples <- numeric(M) # KL para cada muestra MCMC
    for (m in 1:M) {
      kl_samples[m] <- compute_kl(mcmc_array[i, , m], q)
    }
    kl_values[i] <- mean(kl_samples) # Promedio de KL sobre las realizaciones MCMC
  }
  
  return(kl_values) # Retorna los valores KL promedio para cada línea
}



# Function to compute Bhattacharyya distance between two distributions
compute_bhattacharyya <- function(p, q) {
  if (length(p) != length(q)) stop("The distributions p and q must have the same length.")
  
  # Ensure reference distribution is normalized
  q <- q / sum(q)
  
  # Calculate Bhattacharyya coefficient
  bc <- sum(sqrt(p * q))
  
  # Bhattacharyya distance
  d_b <- -log(bc)
  
  return(d_b)
}

# Function to compute the average Bhattacharyya distance for a single trait
compute_bhattacharyya_single_trait <- function(mcmc_array, q) {
  # Normalize the reference distribution q
  q <- q / sum(q)
  
  # Dimensions of the array
  n <- dim(mcmc_array)[1] # Number of lines
  k <- dim(mcmc_array)[2] # Number of categories
  M <- dim(mcmc_array)[3] # Number of MCMC realizations
  
  if (length(q) != k) stop("The length of q must match the number of categories k.")
  
  # Initialize a vector to store the average Bhattacharyya distance for each line
  bhattacharyya_values <- numeric(n)
  
  for (i in 1:n) {
    bhattacharyya_samples <- numeric(M)
    for (m in 1:M) {
      bhattacharyya_samples[m] <- compute_bhattacharyya(mcmc_array[i, , m], q)
    }
    bhattacharyya_values[i] <- mean(bhattacharyya_samples) # Average distance over MCMC samples
  }
  
  return(bhattacharyya_values)
}


# Function to compute the Hellinger distance between two distributions
compute_hellinger <- function(p, q) {
  if (length(p) != length(q)) stop("The distributions p and q must have the same length.")
  p <- sqrt(p)
  q <- sqrt(q)
  h <- sqrt(sum((p - q)^2)) / sqrt(2)
  return(h)
}

# Function to compute the average Hellinger distance for a single trait
compute_hellinger_single_trait <- function(mcmc_array, q) {
  # Normalize the reference distribution q to ensure it sums to 1
  q <- q / sum(q)
  
  # Dimensions of the array
  n <- dim(mcmc_array)[1] # Number of lines
  k <- dim(mcmc_array)[2] # Number of categories
  M <- dim(mcmc_array)[3] # Number of MCMC realizations
  
  if (length(q) != k) stop("The length of q must match the number of categories k.")
  
  # Initialize a vector to store the average Hellinger distance for each line
  hellinger_values <- numeric(n)
  
  for (i in 1:n) {
    hellinger_samples <- numeric(M)
    for (m in 1:M) {
      hellinger_samples[m] <- compute_hellinger(mcmc_array[i, , m], q)
    }
    hellinger_values[i] <- mean(hellinger_samples) # Average distance over MCMC samples
  }
  
  return(hellinger_values)
}


# Main Ordinal Single trait PS function 
OrdinalPS <- function(Xcand, B, thresholds, target = NULL, method = NULL) {
  
  # Revisar que la función de pérdida sea válida
  if (!(method %in% c("kl", "hellinger", "bhattacharyya"))) stop("The method is not valid.\n")
  
  # Threshold matrix
  threshold_matrix <- thresholds
  
  # eta
  eta_matrix <- tcrossprod(Xcand, B)
  
  # Validate dimensions
  if (ncol(Xcand) != ncol(B)) stop("The number of columns in Xcand must match the number of columns in B.")
  if (ncol(eta_matrix) != nrow(thresholds)) stop("The number of columns in thresholds must match the number of rows in eta.")
  
  # Input dimensions
  n_obs <- nrow(eta_matrix)
  n_iter <- ncol(eta_matrix)
  n_thresholds <- ncol(threshold_matrix)
  n_categories <- n_thresholds + 1
  
  if (length(target) != n_categories) stop("The length of target must match the number of categories.")
  
  # Initialize array for probabilities
  probs_array <- array(NA, dim = c(n_obs, n_categories, n_iter))
  
  # Calculate probabilities for each iteration
  for (iter in seq_len(n_iter)) {
    thresholds <- threshold_matrix[iter, ]  # Extract thresholds for the iteration
    for (obs in seq_len(n_obs)) {
      eta <- eta_matrix[obs, iter]  # Extract eta for the observation
      # Calculate cumulative probabilities
      cumulative_probs <- pnorm(thresholds - eta)
      # Calculate category probabilities
      probs_array[obs, 1, iter] <- cumulative_probs[1]
      for (j in 2:n_thresholds) {
        probs_array[obs, j, iter] <- cumulative_probs[j] - cumulative_probs[j - 1]
      }
      probs_array[obs, n_categories, iter] <- 1 - cumulative_probs[n_thresholds]
      # Ensure normalization
      probs_array[obs, , iter] <- probs_array[obs, , iter] / sum(probs_array[obs, , iter])
    }
  }
  
  # Average probabilities across iterations
  yHat <- apply(probs_array, c(1, 2), mean)
  colnames(yHat) <- paste0("Prob_Cat_", 1:ncol(yHat))
  
  # Calculate the distance metric based on the method
  if (method == "kl") {
    e.loss <- compute_kl_single_trait(probs_array, target)
  } else if (method == "hellinger") {
    e.loss <- compute_hellinger_single_trait(probs_array, target)
  } else if (method == "bhattacharyya") {
    e.loss <- compute_bhattacharyya_single_trait(probs_array, target)
  }
  
  ranking <- rank(e.loss)
  
  # Output
  out <- list(method = method, loss = e.loss, ranking = ranking, yHat = data.frame(yHat))
  class(out) <- "MPS"
  
  return(out)
}




