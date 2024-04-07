#' Get average shared alleles with relatives
#'
#' Get average shared alleles of each candidate with the rest of candidates
#'
#' @param Xcand (a matrix) of predictors of the candidates for selection. Dimensions \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set and \eqn{k} is the
#'      number of predictors.
#' @param num_cores (Integer) Number of cores to use.
#' @return A squared matrix of dimension \eqn{n \times n} with proportion of shared alleles with relatives for each candidate of selection.
#' @export 
SimMatrix <- function(Xcand, num_cores = NULL){
  
  `%ni%` <- Negate(`%in%`)
  
  # Input checking
  # Revisamos que no existan valores distintos 0s y 1s.
  if(length(unique(c(Xcand))) %ni% c(2,3) ){
    stop("The matrix of markers should contain only two (dominance) or three (additive) different digits.\n")
  }
  
  n <- nrow(Xcand)
  p <- ncol(Xcand)
  
  # Register parallel backend
  if(is.null(num_cores)) num_cores <- (detectCores() - 1)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  comb <- combn(n, 2)  # Generar todas las comb posibles de filas
  
  # Parallelize the loop using foreach
  out <- matrix(0, n, n)
  
  comb <- foreach(i = 1:ncol(comb), .combine = cbind) %dopar% {
    fila1 <- Xcand[comb[1, i], ]
    fila2 <- Xcand[comb[2, i], ]
    proportion <- mean(fila1 == fila2)
    c(proportion, proportion)  # Fill both upper and lower triangle
  }
  
  v <- comb[1, ]
  out[upper.tri(out, diag = FALSE)] <- v #Upper triangle
  # t(out)   #Lower triangle
  out <- constuyeSimOnes(out)
  
  # Stop the cluster
  stopCluster(cl = cl)
  
  return(out)
}
# SimMatrix <- function(Xcand){
#   
#   `%ni%` <- Negate(`%in%`)
#   
#   # Input checking
#   # Revisamos que no existan valores distintos 0s y 1s.
#   if(length(unique(c(Xcand))) %ni% c(2,3) ) stop("Matrix markers should contain only two (dominance) or three (additive) different digits.\n")
#   
#   n <- nrow(Xcand)
#   p <- ncol(Xcand)
#   
#   comb <- combn(n, 2)  # Generar todas las comb posibles de filas
#   
#   out <- matrix(0, n, n)  # Markers para almacenar los out
#   
#   for (i in 1:ncol(comb)) {
#     fila1 <- Xcand[comb[1, i], ]
#     fila2 <- Xcand[comb[2, i], ]
#     
#     proportion <- mean(fila1 == fila2)  # Comparar las filas y sumar los valores booleanos
#     
#     out[comb[1, i], comb[2, i]] <- proportion
#     out[comb[2, i], comb[1, i]] <- proportion
#     diag(out) <- 1
#   }
#   
#   return(out)
# }
