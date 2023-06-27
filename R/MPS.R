#' Multitrait Parental Selection
#'
#' Function for computing Kullback-Leibler, Energy Score, or Multivariate Asymmetric Loss in Multitrait Parental Selection.
#'
#' @param Xcand (a matrix) of predictors of the candidates for selection. Dimensions \eqn{n \times k}, where
#'      \eqn{n} is the number of individuals in the candidate set and \eqn{k} is the
#'      number of predictors.
#' @param B0 (matrix) of dimension \eqn{M \times t} for the intercept term in the
#'      linear model, where \eqn{M} is the number of markov chain monte carlo samples and \eqn{t} is the number of traits.
#' @param B (array) containing regression coeficients of dimension \eqn{M \times k \times t} of Markov Chain Monte Carlo samples.
#' @param R (matrix) of dimension \eqn{M \times (t \times (t + 1) / 2)} of Markov Chain Monte Carlo samples of the variance-covariance
#'     components in the residual covariance matrix.
#' @param target (vector) of length equal to number of traits (\eqn{t}) reflecting
#'      the breeder's expectation. Default is NULL.
#' @param method (string) the loss function to be used. This must be one of "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".
#' @param measure (string) the distance measure to be used to calculate average distances between lines. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @return A list with Breeding values, posterior expected loss, ranking, and average distance for each candidate of selection.
#' @references 
#'    Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
#'    
#'    Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. 11(2). 1-14
#' @export
#'
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
#' n <- nrow(Y)
#'
#' # Just for example: parental population
#' porc_parental <- 0.4
#' idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
#' XTrn <- X[-idPar,]
#' YTrn <- Y[-idPar,]
#'
#' # ModelFit using BGLR
#' ETA <- list(list(X = XTrn, model = "BRR", saveEffects = TRUE))
#' model <- Multitrait(y = YTrn,
#'                   ETA = ETA,
#'                   intercept = TRUE,
#'                   resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = TRUE),
#'                   saveAt = paste0("./chains/", "BRR_"),
#'                   nIter = 50000,
#'                   burnIn = 20000)
#'
#' # Reading Posterior MCMC
#' mu <- as.matrix(read.table(file = "./chains/BRR_mu.dat", header = FALSE))
#' B <- readBinMatMultitrait('./chains/BRR_ETA_1_beta.bin')
#' R <- as.matrix(read.table(file = "./chains/BRR_R.dat", header = FALSE))
#' XPar <- X[idPar, ]
#'
#' # Evaluating Loss Function
#' out <- MPS(Xcand = XPar,
#'            B0 = mu,
#'            B = B,
#'            R = R,
#'            p = 0.1,
#'            method = "kl")
#'
#' # Plotting results
#' colnames(out$yHat) <- colnames(Y)
#' pairs(out$yHat,
#'      col = ifelse(out$selected, "red", "darkgray"),
#'      pch = ifelse(out$selected, 19, 1))
#' }
EvalMPS <- function(Xcand, B0, B, R, target = NULL, method = "kl", measure = "euclidean", verbose = FALSE){

  if(!is.matrix(Xcand)) stop("Xcand must be a matrix.\n")
  if(!is.matrix(B0)) stop("B0 must be a matrix.\n")
  if(!is.array(B)) stop("B must be an array.\n")
  if(!is.matrix(R)) stop("R must be a matrix.\n")
  # if(p <= 0 | p >= 1) stop("p should be a number between 0 and 1.\n")
  if(!(method %in% c("kl", "malf", "energy"))) stop("The method is not valid.\n")
  if(!(measure %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))) stop("The method is not valid.\n")

  t <- ncol(B0)
  if(t < 2) stop("The number of traits must be at least two.\n")

  # Start calculations
  means <- colMeans(B0)
  desvEst <- sqrt(diag(xpnd(colMeans(R))))
  
  if(is.null(target)){
    target <- means + 2*desvEst  
  }else{
    if(!is.numeric(target)) stop("target must be a real vector.\n")
    if(t != length(target)) stop("The length of target must be the same the number of traits.\n")
  }

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

  # Average distances
  aveDistances <- aveDist(Xcand, measure)

  # Calculating posterior expected loss and Breeding values
  # nSelected <- ceiling(n*p)
  e.loss <- apply(loss, 1, mean, na.rm = TRUE)
  ranking <- rank(e.loss)
  # selected <- order(e.loss, decreasing = FALSE)[1:nSelected]
  # selected <- ifelse(1:n %in% selected, TRUE, FALSE)
  yHat <- apply(mu_c, c(1,2), mean)
  out <- list(method = method, loss = e.loss, ranking = ranking, aveDist = aveDistances, yHat = data.frame(yHat))
  class(out) <- "MPS"
  return(out)
}
