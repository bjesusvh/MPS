#' Approximation in Multitrait Parental Selection
#'
#' Approximate the Kullback Leibler, Energy Score or Multivariate Assymetric Loss.
#'
#' @param B0 (vector) with the overall mean of each trait of length equal to number of traits (\eqn{t}).
#' @param yHat (matrix) of posterior mean of each candidate of dimension \eqn{n \times t}, where
#'      \eqn{n} is the number of individuals in the candidate set.
#' @param R (square matrix) of the residual covariance matrix of dimension \eqn{t \times t}.
#' @param target (vector) of length equal to number of traits (\eqn{t}) reflecting
#'      the breeder's expectation. Default is NULL.
#' @param method (string) the loss function to be used. This must be one of "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".
#' @param direction A vector of length equal to number of traits (\eqn{t}) reflecting the direction
#'     of improvement desired, 1 for increase goal, -1 for decreasing goal.
#'     Default is 1 (increase) for all traits.
#' @return A list with Breeding values, expected loss, and the rank for each candidate of selection.
#' @references 
#'    Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
#'    
#'    Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. 11(2). 1-14
#' @export
#' @examples
#' \dontrun{
#'
#' }
ApproxMPS <- function(B0, yHat, R, target = NULL, method = "kl", direction = NULL){

  # Checking inputs
  if(!is.vector(B0)) stop("B0 must be a vector.\n")
  if(!is.matrix(yHat)) stop("yHat must be a matrix.\n")
  if(!is.matrix(R)) stop("R must be a matrix.\n")
  if(!(method %in% c("kl", "malf", "energy"))) stop("The method is not valid.\n")

  t <- length(B0)
  if(t < 2) stop("The number of traits must be at least two.\n")
  if(is.null(direction)){
    direction <- rep(1, t)
  }else{
    if(t != length(direction)) stop("The length of direction must be the same the number of traits.\n")
    if(!(sum(direction %in% c(-1,1)) == t)) stop("Direction should be -1 or 1.\n")
  }

  # Start calculations
  desvEst <- sqrt(diag(R))

  # Target
  if(is.null(target)){
    target <- B0 + (2*desvEst*direction)
  }else{
    if(!is.numeric(target)) stop("target must be a real vector.\n")
    if(t != length(target)) stop("The length of target must be the same the number of traits.\n")
  }

  # mvtnorm::
  muS <- as.vector(mtmvnorm(lower = target,
                            upper = rep(Inf, t),
                            mean = B0,
                            sigma = R,
                            doComputeVariance = FALSE)$tmean)

  if(method == "kl"){
    Kinv <- solve(R)
    e.loss <- aproxMultiKL(target, rep(Inf, t), B0, yHat, muS, R, Kinv)
  }else if(method == "energy"){
    e.loss <- aproxEnergyScore(yHat, muS)
  }
  else{
    tau <- rep(0.95, t)
    e.loss <- aproxMALF(yHat, muS, tau)
  }
  
  # Calculating posterior expected loss and Breeding values
  # n <- nrow(yHat)
  # nSelected <- ceiling(n*p)
  ranking <- rank(e.loss)
  # selected <- order(e.loss, decreasing = FALSE)[1:nSelected]
  # selected = ifelse(1:n %in% selected, TRUE, FALSE)
  out <- list(method = method, loss = e.loss, ranking = ranking, yHat = data.frame(yHat))
  class(out) <- "MPS"
  return(out)
}



