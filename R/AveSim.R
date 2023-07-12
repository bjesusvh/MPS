#' Get average similarities
#'
#' Get average similarities for each candidate of selection
#'
#' @param K The pedigree, genomic relationship or similarity matrix. Dimensions \eqn{n \times n}, where
#'      \eqn{n} is the number of individuals in the candidate set.
#' @return A vector with average similarities for each candidate of selection.
#' @export 
#'
#' @examples
#' \dontrun{
#' set.seed(624) # for reproducibility
#' X <- matrix(sample(c(-1,0,1), 8*4, replace = TRUE), nrow = 4)
#' simM <- SimMatrix(X)
#' AveSim(simM)
#' }
AveSim <- function(K){ 
  
  if(dim(K)[1] != dim(K)[2]) stop("K should be a square matrix.")
  
  # Computations
  diag(K) <- NA
  out <- apply(K, 1, function(x) mean(x, na.rm = TRUE))
  
  # return
  return(out)
}


# getSim <- function(K){ 
#   
#   if(dim(K)[1] != dim(K)[2]) stop("K should be a square matrix.")
#   
#   # Computations
#   upper_tri <- upper.tri(K)
#   indices <- which(upper_tri, arr.ind = TRUE)
#   values <- K[indices]
#   df <- data.frame(row = indices[, 1], 
#                    col = indices[, 2],
#                    value = values)
#   # Mean values
#   rel1 <- tapply(df$value, INDEX = factor(df$row), FUN = mean)
#   rel2 <- tapply(df$value, INDEX = factor(df$col), FUN = mean)
#   aveSim <- c(rel1, rel2[length(rel2)])
#   return(aveSim)
# }