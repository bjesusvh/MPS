% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AveSim.R
\name{AveSim}
\alias{AveSim}
\title{Get average similarities}
\usage{
AveSim(K)
}
\arguments{
\item{K}{The pedigree, genomic relationship or similarity matrix. Dimensions \eqn{n \times n}, where
\eqn{n} is the number of individuals in the candidate set.}
}
\value{
A vector with average similarities for each candidate of selection.
}
\description{
Get average similarities for each candidate of selection
}
\examples{
\dontrun{
set.seed(624) # for reproducibility
X <- matrix(sample(c(-1,0,1), 8*4, replace = TRUE), nrow = 4)
simM <- SimMatrix(X)
AveSim(simM)
}
}
