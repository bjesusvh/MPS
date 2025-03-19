#' Single-Trait Ordinal Parental Selection
#'
#' Function to compute Kullback-Leibler, Hellinger, or Bhattacharyya loss in Single-Trait Ordinal Parental Selection.
#'
#' @param Xcand A matrix of predictors for the candidates under selection, with dimensions \eqn{n \times k}.
#'    \eqn{n}: Number of candidates.
#'    \eqn{k}: Number of predictors.
#' @param B A matrix containing regression coefficients, with dimensions \eqn{M \times k}.
#'    \eqn{M}: Number of Markov Chain Monte Carlo (MCMC) samples.
#'    \eqn{k}: Number of predictors.
#' @param thresholds A matrix of thresholds values from ordinal regression, with dimensions \eqn{M \times t}, derived from Markov Chain Monte Carlo (MCMC) samples.
#'    Thresholds are the values that define the boundaries between the ordinal categories in the outcome variable. 
#'    They are cut-points on the underlying scale of the latent variable that determine the category boundaries
#' @param target A vector that contain probabilities or proportions for the categories within that trait. 
#'    The probabilities/proportions must sum to 1. Zero entries are not allowed.
#' @param method A string indicating the loss function to be used. Options are "kl" for Kullback-Leibler, "hellinger" for Hellinger, or "bhattacharyya". Default is "kl".
#' @return A list containing the posterior expected loss, the rank for each candidate, and the probability of each category.
#' @references 
#'    Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
#'    
#'    Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. G3 Genes|Genomes|Genetics. 11(2). 1-14
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of OrdinalMPS
#' }
OrdinalPS <- function(Xcand, B, thresholds, target = NULL, method = "kl") {
  
  # Revisar que la función de pérdida sea válida
  if (!(method %in% c("kl", "hellinger", "bhattacharyya"))) stop("The method is not valid.\n")
  
  # Matriz de umbrales
  threshold_matrix <- thresholds
  
  # eta
  eta_matrix <- tcrossprod(Xcand, B)
  
  # Validar dimensiones
  if (ncol(Xcand) != ncol(B)) stop("The number of columns in Xcand must match the number of columns in B.")
  if (ncol(eta_matrix) != nrow(thresholds)) stop("The number of columns in thresholds must match the number of rows in eta.")
  
  # Dimensiones de entrada
  n_obs <- nrow(eta_matrix)
  n_iter <- ncol(eta_matrix)
  n_thresholds <- ncol(threshold_matrix)
  n_categories <- n_thresholds + 1
  
  if (length(target) != n_categories) stop("The length of target must match the number of categories.")
  
  # Inicializa el arreglo de probabilidades
  probs_array <- array(NA, dim = c(n_obs, n_categories, n_iter))
  
  # Calcular las probabilidades para cada iteración
  for (iter in seq_len(n_iter)) {
    thresholds <- threshold_matrix[iter, ]  # Extrae los umbrales para la iteración
    for (obs in seq_len(n_obs)) {
      eta <- eta_matrix[obs, iter]  # Extrae eta para la observación
      # Calcular probabilidades acumuladas
      cumulative_probs <- pnorm(thresholds - eta)
      # Calcular probabilidades por categoría
      probs_array[obs, 1, iter] <- cumulative_probs[1]
      for (j in 2:n_thresholds) {
        probs_array[obs, j, iter] <- cumulative_probs[j] - cumulative_probs[j - 1]
      }
      probs_array[obs, n_categories, iter] <- 1 - cumulative_probs[n_thresholds]
      # Asegurar normalización
      probs_array[obs, , iter] <- probs_array[obs, , iter] / sum(probs_array[obs, , iter])
    }
  }
  
  # Promedio de probabilidades sobre iteraciones
  yHat <- apply(probs_array, c(1, 2), mean)
  colnames(yHat) <- paste0("Prob_Cat_", 1:ncol(yHat))
  #result <- apply(yHat, 1, which.max)
  
  # Calcular la métrica de distancia según el método
  e.loss <- compute_divergence_single_trait(probs_array, target, method)
  
  ranking <- rank(e.loss)
  
  # Salida
  out <- list(method = method, loss = e.loss, ranking = ranking, yHat = data.frame(yHat))
  class(out) <- "MPS"
  
  return(out)
}