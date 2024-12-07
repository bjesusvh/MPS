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
