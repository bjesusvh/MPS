% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultitraitOPS.R
\name{MultitraitOPS}
\alias{MultitraitOPS}
\title{Multitrait Ordinal Parental Selection}
\usage{
MultitraitOPS(Xcand, B, thresholds, target = NULL, method = "kl")
}
\arguments{
\item{Xcand}{A matrix of predictors for the candidates under selection, with dimensions \eqn{n \times k}.
\eqn{n}: Number of candidates.
\eqn{k}: Number of predictors.}

\item{B}{A list of length equal to the number of traits. Each element in the list is a matrix containing regression coefficients, with dimensions \eqn{M \times k}.
\eqn{M}: Number of Markov Chain Monte Carlo (MCMC) samples.
\eqn{k}: Number of predictors.}

\item{thresholds}{A list of matrices that contains thresholds values from ordinal regression, with dimensions \eqn{M \times t}, derived from Markov Chain Monte Carlo (MCMC) samples.
Thresholds are the values that define the boundaries between the ordinal categories in the outcome variable. 
They are cut-points on the underlying scale of the latent variable that determine the category boundaries}

\item{target}{A list of vectors, where each vector corresponds to an ordinal trait. Each vector must contain probabilities or proportions for the categories within that trait. 
The length of each vector must match the number of categories in the trait, and the probabilities/proportions must sum to 1. Zero entries are not allowed.}

\item{method}{A string indicating the loss function to be used. Options are "kl" for Kullback-Leibler, "hellinger" and "bhattacharyya". Default is "kl".}
}
\value{
A list containing the expected loss (distance), the rank for each candidate, and the predicted category for each trait.
}
\description{
Function to compute Kullback-Leibler, Hellinger, or Bhattacharyya loss in Multitrait Ordinal Parental Selection.
}
\examples{
\dontrun{
# Example usage of OrdinalMPS
}
}
\references{
Villar-Hernández, B.J., et.al. (2018). A Bayesian Decision Theory Approach for Genomic Selection. G3 Genes|Genomes|Genetics. 8(9). 3019–3037
   
   Villar-Hernández, B.J., et.al. (2021). Application of multi-trait Bayesian decision theory for parental genomic selection. G3 Genes|Genomes|Genetics. 11(2). 1-14
}
