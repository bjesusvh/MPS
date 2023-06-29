#' Approximation in Multitrait Parental Selection
#'
#' Approximate the Kullback Leibler, Energy Score or Multivariate Assymetric Loss.
#'
#' @param B0 (vector) with the overall mean of each trait of length equal to number of traits (\eqn{t}).
#' @param ETA (matrix) of posterior mean of each candidate of dimension \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set.
#' @param R (square matrix) of the residual covariance matrix of dimension \eqn{t \times t}.
#' @param target (vector) of length equal to number of traits (\eqn{t}) reflecting
#'      the breeder's expectation. Default is NULL.
#' @param method (string) the loss function to be used. This must be one of "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".
#' @param A (matrix) pedigree information of dimension \eqn{n \times n}.
#' @return A list with Breeding values, posterior expected loss, and ranking for each candidate of selection.
#' @references 
#'    Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
#'    
#'    Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. 11(2). 1-14
#' @export
#' @examples
#' \dontrun{
#' # Clean and setting local environment
#' rm(list = ls())
#'
#' setwd("Put your working directory here")
#'
#' # Loading needed packages
#' library(MPS)
#' library(BGLR)
#'
#' # Loading dataset
#' data(wheat)
#' Y <- as.matrix(wheat.Y)
#' X <- scale(wheat.X, center=TRUE, scale = TRUE)
#' K <- wheat.A
#' n <- nrow(Y)
#'
#' # Just for example: parental population
#' idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
#' YTrn <- Y
#' YTrn[idPar,] <- NA
#'
#' # ModelFit using BGLR
#' ETA <- list(list(X = X, model = "BRR"),
#'             list(K = K, model = "RKHS"))
#'
#' model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,
#'                     resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = FALSE),
#'                     nIter = 30000,
#'                     burnIn = 5000)
#'
#' # Retrieving puntual estimates
#' R <- as.matrix(model$resCov$R)               # residual cov matrix
#' B0 <- as.numeric(model$mu)                   # Overall mean
#' yHat <- model$ETAHat[model$missing_records,] # Puntual BVs
#'
#' # Eval loss functions
#' out <- ApproxMPS(B0 = B0, ETA = yHat, R = R, p = 0.1, method = "kl")
#'
#' # Plotting results
#' pairs(out$yHat,
#'   col = ifelse(out$selected, "red", "darkgray"),
#'   pch = ifelse(out$selected, 19, 1))
#' }
ApproxMPS <- function(B0, ETA, R, target = NULL, method = "kl", A = NULL){

  # Checking inputs
  if(!is.vector(B0)) stop("B0 must be a vector.\n")
  if(!is.matrix(ETA)) stop("ETA must be a matrix.\n")
  if(!is.matrix(R)) stop("R must be a matrix.\n")
  # if(p <= 0 | p >= 1) stop("p should be a number between 0 and 1.\n")
  if(!(method %in% c("kl", "malf", "energy"))) stop("The method is not valid.\n")

  t <- length(B0)
  if(t < 2) stop("The number of traits must be at least two.\n")

  # Start calculations
  desvEst <- sqrt(diag(R))

  # Target
  if(is.null(target)){
    target <- B0 + 2*desvEst
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
    e.loss <- aproxMultiKL(target, rep(Inf, t), B0, ETA, muS, R, Kinv)
  }else if(method == "energy"){
    e.loss <- aproxEnergyScore(ETA, muS)
  }
  else{
    tau <- rep(0.95, t)
    e.loss <- aproxMALF(ETA, muS, tau)
  }
  
  # Calculating similarities
  if(!is.null(A)){
    aveSimilarities <- aveSim(A)
  }else{
    aveSimilarities <- NULL
  }

  # Calculating posterior expected loss and Breeding values
  # n <- nrow(ETA)
  # nSelected <- ceiling(n*p)
  ranking <- rank(e.loss)
  # selected <- order(e.loss, decreasing = FALSE)[1:nSelected]
  # selected = ifelse(1:n %in% selected, TRUE, FALSE)
  out <- list(method = method, loss = e.loss, ranking = ranking, aveSim = aveSimilarities, yHat = data.frame(ETA))
  class(out) <- "MPS"
  return(out)
}



