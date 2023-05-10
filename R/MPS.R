#' Multitrait Parental Selection
#'
#' multiLoss function to compute Kullback Leibler, Energy Score or Multivariate Assymetric Loss in Multitrait Parental Selection.
#'
#' @param Xcand (a matrix) of predictors of candidates for selection of dimensions \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set and \eqn{k} is the
#'      number of predictors.
#' @param B0 (matrix) of dimension \eqn{M \times t} for the intercept term in
#'      linear model, where \eqn{M} is the number of markov chain monte carlo samples and \eqn{t} is the number of traits.
#' @param B (array) contain regression coeficients of dimension, \eqn{M \times k \times t} of Markov Chain Monte Carlo samples.
#' @param R (matrix) of dimension \eqn{M \times (t \times (t + 1) / 2)} of Markov Chain Monte Carlo samples of the variance-covariance
#'     components in the residual covariance matrix.
#' @param p (escalar) the proportion of selected candidates.
#' @param method (string) indicate the loss function to use. Posible values are "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".
#' @param verbose (logical) if TRUE the iteration history is printed, default FALSE.
#' @return A list of BVs, posterior expected loss and a logical vector that indicate what lines are selected.
#' @export
#'
#' @examples
#' \dontrun{
#' # Clean and setting local environment
#' rm(list = ls())
#' setwd("Put your working directory here")
#' # Loading needed packages
#' library(scoringGS)
#' library(BGLR)
#' # Loading dataset
#' data(wheat)
#' Y <- as.matrix(wheat.Y)
#' X <- scale(wheat.X, center=TRUE, scale = TRUE)
#' n <- nrow(Y)
#' # Just for example: parental population
#' porc_parental <- 0.4
#' idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
#' XTrn <- X[-idPar,]
#' YTrn <- Y[-idPar,]
#' # ModelFit using BGLR
#' ETA <- list(list(X = XTrn, model = "BRR", saveEffects = TRUE))
#' model <- Multitrait(y = YTrn,
#'                   ETA = ETA,
#'                   intercept = TRUE,
#'                   resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = TRUE),
#'                   saveAt = paste0("./chains/", "BRR_"),
#'                   nIter = 50000,
#'                  burnIn = 20000)
#' # Reading Posterior MCMC
#' mu <- as.matrix(read.table(file = "./chains/BRR_mu.dat", header = FALSE))
#' B <- readBinMatMultitrait('./chains/BRR_ETA_1_beta.bin')
#' R <- as.matrix(read.table(file = "./chains/BRR_R.dat", header = FALSE))
#' XPar <- X[idPar, ]
#' # Evaluating Loss Function
#' out <- MPS(Xcand = XPar,
#'            B0 = mu,
#'            B = B,
#'            R = R,
#'            p = 0.1,
#'            method = "kl")
#' # Plotting results
#' colnames(out$yHat) <- colnames(Y)
#' pairs(out$yHat,
#'      col = ifelse(out$selected, "#2c7fb8", "darkgray"),
#'      pch = ifelse(out$selected, 19, 1))
#' }
MPS <- function(Xcand, B0, B, R, p, method = "kl", verbose = FALSE){

  #Xcand, B0, B, R, target, p, method = "kl", verbose = FALSE
  #@param target A vector of length equal to number of traits (\eqn{t}) reflecting
  #     the breeder's expectation, i.e., the gain they expect to observe after selection.
  # Checking inputs
  if(!is.matrix(Xcand)) stop("Xcand must be a matrix.\n")
  if(!is.matrix(B0)) stop("B0 must be a matrix.\n")
  if(!is.array(B)) stop("B must be an array.\n")
  if(!is.matrix(R)) stop("R must be a matrix.\n")
  # if(!is.numeric(target)) stop("target must be a real vector.\n")
  if(p <= 0 | p >= 1) stop("p should be a number between 0 and 1.\n")
  if(!(method %in% c("kl", "malf", "energy"))) stop("The method is not valid.\n")

  t <- ncol(B0)
  if(t < 2) stop("The number of traits must be at least two.\n")
  # if(t != length(target)) stop("The length of target must be the same the number of traits.\n")

  # Start calculations
  means <- colMeans(B0)
  desvEst <- sqrt(diag(xpnd(colMeans(R))))
  target <- means + 2*desvEst
  Pall <- lapply(1:nrow(R), FUN = toNearPD, x = R)

  if(method == "malf"){
    tau <- rep(0.95, t)
  }

  if(method == "kl"){
    Pallinv <- lapply(Pall, solve)
  }

  mu_c <- posteriorBV(Xcand, B, B0)
  n <- dim(mu_c)[1]
  M <- dim(mu_c)[3]
  loss <- matrix(NA, n, M)

  # Evaluation of loss function using MCMC

  for(i in 1:M) {

    time <- proc.time()[3]
    mu1 <- as.vector(B0[i,])
    mu2 <- mu_c[,,3]
    upper <- rep(Inf, length(target))
    K <- Pall[[i]]

    # mvtnorm::
    muS <- as.vector(mtmvnorm(lower = target,
                              upper = upper,
                              mean = mu1,
                              sigma = K,
                              doComputeVariance = FALSE)$tmean)

    # Choose the loss function to evaluate
    if(method == "kl"){
      Kinv <- Pallinv[[i]]
      loss[,i] <- multi_KL(target, upper, mu1, mu2, muS, K, Kinv)
    }else if(method == "energy"){
      loss[,i] <- EnergyScore(mu2, muS, K)
    }
    else{
      loss[,i] <- MALF(mu2, muS, K, tau)
    }

    if(verbose){
      tmp <- proc.time()[3]
      cat(c(paste(c("  Iter = ", "Time/Iter = "), round(c(i, c(tmp - time)), 2), sep = ""), "s"), "\n")
      time <- tmp
    }
  }

  # Calculating posterior expected loss and Breeding values
  nSelected <- ceiling(n*p)
  e.loss <- apply(loss, 1, mean, na.rm = TRUE)
  selected <- order(e.loss, decreasing = FALSE)[1:nSelected]
  yHat <- apply(mu_c, c(1,2), mean)
  out <- list(method = method, loss = e.loss, yHat = data.frame(yHat), selected = ifelse(1:n %in% selected, TRUE, FALSE))
  class(out) <- "MPS"
  return(out)
}
