#' Multitrait Ordinal Parental Selection
#'
#' Function to compute Kullback-Leibler, Hellinger, or Bhattacharyya loss in Multitrait Ordinal Parental Selection.
#'
#' @param Xcand A matrix of predictors for the candidates under selection, with dimensions \eqn{n \times k}.
#'    \eqn{n}: Number of candidates.
#'    \eqn{k}: Number of predictors.
#' @param B A list of length equal to the number of traits. Each element in the list is a matrix containing regression coefficients, with dimensions \eqn{M \times k}.
#' - \eqn{M}: Number of Markov Chain Monte Carlo (MCMC) samples.
#' - \eqn{k}: Number of predictors.
#' @param thresholds A list of matrices that contains thresholds values from ordinal regression, with dimensions \eqn{M \times t}, derived from Markov Chain Monte Carlo (MCMC) samples.
#'    Thresholds are the values that define the boundaries between the ordinal categories in the outcome variable. 
#'    They are cut-points on the underlying scale of the latent variable that determine the category boundaries
#' @param target A list of vectors, where each vector corresponds to an ordinal trait. Each vector must contain probabilities or proportions for the categories within that trait. 
#' The length of each vector must match the number of categories in the trait, and the probabilities/proportions must sum to 1. Zero entries are not allowed.
#' @param method A string indicating the loss function to be used. Options are "kl" for Kullback-Leibler, "hellinger" for Hellinger, or "bhattacharyya". Default is "kl".
#' @return A list containing the posterior expected loss, the rank for each candidate, and the predicted category for each trait.
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
MultitraitOPS <- function(Xcand, B, thresholds, target = NULL, method = "kl") {
  
  # Validación del método
  if (!(method %in% c("kl", "hellinger", "bhattacharyya"))) {
    stop("Invalid method. Choose 'kl', 'hellinger', or 'bhattacharyya'.")
  }

  # Validación de entradas no nulas
  if (is.null(target)) stop("Argument 'target' must be a non-null list of vectors.")
  if (!is.matrix(Xcand)) stop("Argument 'Xcand' must be a matrix.")
  if (!is.list(B)) stop("Argument 'B' must be a list of matrices.")
  if (!is.list(thresholds)) stop("Argument 'thresholds' must be a list of matrices.")
  if (!is.list(target)) stop("Argument 'target' must be a list of vectors.")

  # Validación de dimensiones de traits
  Ntraits_B <- length(B)
  Ntraits_th <- length(thresholds)
  Ntraits_target <- length(target)
  
  if (!(Ntraits_B == Ntraits_th && Ntraits_th == Ntraits_target)) {
    stop("Mismatch in the number of traits among 'B', 'thresholds', and 'target'.")
  }
  
  # Procesamiento por cada trait
  out <- lapply(seq_len(Ntraits_B), function(t) {
    OrdinalPS(Xcand, B[[t]], thresholds[[t]], target[[t]], method = method)
  })
  
  names(out) <- paste0("Trait ", seq_len(Ntraits_B))
  
  # Cálculo de pérdida total
  e.loss <- Reduce(`+`, lapply(out, `[[`, 2))
  
  # Predicción de categorías para cada trait
  result <- lapply(out, function(trait) apply(trait$yHat, 1, which.max))
  
  # Ranking global
  ranking <- rank(e.loss)
  
  # Salida
  salida <- list(loss = e.loss, ranking = ranking, Cat_Predicted = result)
  
  return(salida)
}