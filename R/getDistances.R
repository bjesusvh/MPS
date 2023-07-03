#' Get average distances
#'
#' Get average distance for each candidate of selection
#'
#' @param Xcand (a matrix) of predictors of the candidates for selection. Dimensions \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set and \eqn{k} is the
#'      number of predictors.
#' @param measure (string) the distance measure to be used to calculate average distances between lines. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @return A vector with average distance for each candidate of selection.
#' @export 
#'
#' @examples
#' \dontrun{
#' X <- matrix(rbinom(25, 1, p=0.4), nrow = 5, ncol = 5)
#' getDist(X)
#' }
getDist <- function(Xcand, measure = "euclidean"){ 
  
  if(!(measure %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))) stop("The method is not valid.\n")
  
  D <- as.matrix(dist(Xcand, method = measure))
  # Computations
  upper_tri <- upper.tri(D)
  indices <- which(upper_tri, arr.ind = TRUE)
  values <- D[indices]
  df <- data.frame(row = indices[, 1], 
                   col = indices[, 2],
                   value = values)
  # Mean values
  rel1 <- tapply(df$value, INDEX = factor(df$row), FUN = mean)
  rel2 <- tapply(df$value, INDEX = factor(df$col), FUN = mean)
  aveDist <- c(rel1, rel2[length(rel2)])
  return(aveDist)
}