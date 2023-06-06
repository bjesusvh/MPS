#' Multitrait Parental Selection
#'
#' Approximate the Kullback Leibler, Energy Score or Multivariate Assymetric Loss.
#'
#' @param B0 (vector) with the overall mean of each trait of length equal to number of traits (\eqn{t}).
#' @param ETA (matrix) of posterior mean of each candidate of dimension \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set.
#' @param R (square matrix) of the residual covariance matrix of dimension \eqn{t \times t}.
#' @param p (escalar) the proportion of selected candidates.
#' @param method (string) indicate the loss function to use. Posible values are "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".
#' @return A list of BVs, approximated expected loss and a logical vector that indicate what lines are selected.
#' @references 
#'    Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
#' @export
#' @examples
#' \dontrun{
#' # Clean and setting local environment
#' rm(list = ls())
#' setwd("Put your working directory here")
#' # Loading needed packages
#' library(MPS)
#' library(BGLR)
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
ApproxMPS <- function(B0, ETA, R, p, method = "kl"){

  # Checking inputs
  if(!is.vector(B0)) stop("B0 must be a vector.\n")
  if(!is.matrix(ETA)) stop("ETA must be a matrix.\n")
  if(!is.matrix(R)) stop("R must be a matrix.\n")
  # if(!is.numeric(target)) stop("target must be a real vector.\n")
  if(p <= 0 | p >= 1) stop("p should be a number between 0 and 1.\n")
  if(!(method %in% c("kl", "malf", "energy"))) stop("The method is not valid.\n")

  t <- length(B0)
  if(t < 2) stop("The number of traits must be at least two.\n")
  # if(t != length(target)) stop("The length of target must be the same the number of traits.\n")

  # Start calculations
  desvEst <- sqrt(diag(R))
  target <- B0 + 2*desvEst

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

  # Calculating posterior expected loss and Breeding values
  n <- nrow(ETA)
  nSelected <- ceiling(n*p)
  selected <- order(e.loss, decreasing = FALSE)[1:nSelected]
  out <- list(method = method, loss = e.loss, yHat = data.frame(ETA), selected = ifelse(1:n %in% selected, TRUE, FALSE))
  class(out) <- "MPS"
  return(out)
}



