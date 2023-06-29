# MPS
An R Package to Assist Multi-trait Parental Selection via Genomic, Pedigree, or Phenotypic Selection.

## Installation from GitHub:

```
#1 install devtools
install.packages(pkg='devtools',repos='https://cran.r-project.org/')

#2 load the library
library(devtools)

#3 install MPS from GitHub                                          
install_github('https://github.com/bjesusvh/MPS')                     
```

## Examples

**Example using Genomic Selection**

```r
# Cleaning and setting local environment
rm(list = ls())
setwd("path")
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
out <- FastMPS(Xcand = XPar, B0 = B0, B = B, R = R, method = "kl")
```


**Example using Genomic + Pedigree**

**Example using Phenotypic selection**

**Example for Multi-environment best parent selection**
