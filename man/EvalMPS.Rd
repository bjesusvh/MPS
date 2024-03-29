% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvalMPS.R
\name{EvalMPS}
\alias{EvalMPS}
\title{Multitrait Parental Selection}
\usage{
EvalMPS(
  Xcand,
  B0,
  B,
  R,
  target = NULL,
  method = "kl",
  direction = NULL,
  verbose = FALSE
)
}
\arguments{
\item{Xcand}{(a matrix) of predictors of the candidates for selection. Dimensions \eqn{n \times k}, where
\eqn{n} is the number of individuals in the candidate set and \eqn{k} is the
number of predictors.}

\item{B0}{(matrix) of dimension \eqn{M \times t} for the intercept term in the
linear model, where \eqn{M} is the number of markov chain monte carlo samples and \eqn{t} is the number of traits.}

\item{B}{(array) containing regression coefficients of dimension \eqn{M \times k \times t} of Markov Chain Monte Carlo samples.}

\item{R}{(matrix) of dimension \eqn{M \times (t \times (t + 1) / 2)} of Markov Chain Monte Carlo samples of the variance-covariance
components in the residual covariance matrix.}

\item{target}{(vector) of length equal to number of traits (\eqn{t}) reflecting
the breeder's expectation. Default is NULL.}

\item{method}{(string) the loss function to be used. This must be one of "kl" for Kullback-Leibler, "energy" for Energy Score and "malf" for Multivariate Assymetric Loss. Default is "kl".}

\item{direction}{A vector of length equal to number of traits (\eqn{t}) reflecting the direction
of improvement desired, 1 for increase goal, -1 for decreasing goal.
Default is 1 (increase) for all traits.}

\item{verbose}{(logical) if TRUE the iteration history is printed, default FALSE}
}
\value{
A list with Breeding values, posterior expected loss, and the rank for each candidate of selection.
}
\description{
Function for computing Kullback-Leibler, Energy Score, or Multivariate Asymmetric Loss in Multitrait Parental Selection.
}
\examples{
\dontrun{
# Cleaning and setting local environment
rm(list = ls())
setwd("Put your working directory")
library(MPS)
library(BGLR)
# Loading dataset
data(wheat)
Y <- as.matrix(wheat.Y)
X <- scale(wheat.X, center=TRUE, scale = TRUE)
n <- nrow(Y)
# Just for example
set.seed(647) # for reproducibility
porc_parental <- 0.4
idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
XTrn <- X[-idPar,]
YTrn <- Y[-idPar,]
XPar <- X[idPar, ]

# Model fitting using BGLR
ETA <- list(list(X = XTrn, model = "BRR", saveEffects = TRUE))
model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,
                    resCov = list(type = "UN", saveEffects = TRUE), nIter = 100000, burnIn = 30000)

# Reading Posterior MCMC
B0 <- as.matrix(read.table(file = "mu.dat", header = FALSE))
B <- readBinMatMultitrait('ETA_1_beta.bin')
R <- as.matrix(read.table(file = "R.dat", header = FALSE))

# Evaluating Loss Function
out <- EvalMPS(Xcand = XPar, B0 = B0, B = B, R = R, method = "kl")

# Get average distances for each candidate
d <- getDist(XPar)

## Plotting similarity
library(utilidades) # install_github('bjesusvh/utilidades')
plotMPS(MPSObject = out, Dist = d, qLoss = 0.25, qDistance = 0.25)
}
}
\references{
Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
   
   Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. 11(2). 1-14
}
