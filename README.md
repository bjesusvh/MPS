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

**Example 1: Multi-trait Genomic selection**

```r
rm(list = ls())                  # Cleaning R memory
setwd("path in your computer")   # Working directory  
library(MPS); library(BGLR)      # Needed R Packages
data(wheat)                      # Loading dataset
Y <- as.matrix(wheat.Y)          # Phenotypic records  
X <- scale(wheat.X, center=TRUE, scale = TRUE) # Genotypic data
set.seed(647)             # Seed for reproducibility
n <- nrow(Y)              # Number of records
porc_parental <- 0.4      # Percentaje in candidates
idPar <- sort(sample(1:n, # Random observation for candidate set
                     ceiling(porc_parental*n), 
                     replace = FALSE))
XTrn <- X[-idPar,]        # Genotypic data in Training set
XPar <- X[idPar, ]        # Genotypic data in Parental set
YTrn <- Y[-idPar,]        # Phenotypic data in Training set
# Specification of the linear predictor for BGLR
ETA <- list(list(X = X, model = "BRR", saveEffects = TRUE))
# Run the model using Multitrait function of BGLR R Package
model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,   
                    resCov = list(type = "UN", saveEffects = TRUE), 
                    nIter = 100000, burnIn = 30000)
# Reading Posterior MCMC of model parameters to compute expected loss
B0 <- as.matrix(read.table(file = "mu.dat",  # Overall mean
                           header = FALSE))
B <- readBinMatMultitrait('ETA_1_beta.bin')  # Regression coefficients
R <- as.matrix(read.table(file = "R.dat",    # Residual covariance matrix
                          header = FALSE))
# Evaluate loss function
out <- FastMPS(Xcand = XPar, B0 = B0, B = B, R = R, method = "kl")
# Get average distances for each candidate
s <- AveSim(SimMatrix(wheat.X[idPar,]))
out
```

**Example 2: Multi-trait Genomic + pedigree selection**

```r
# Cleaning and setting local environment
rm(list = ls())
setwd("your path")
library(MPS); library(BGLR)

# Loading dataset
data(wheat); Y <- as.matrix(wheat.Y)
X <- scale(wheat.X, center=TRUE, scale = TRUE)
K <- wheat.A; n <- nrow(Y)

# Just for example  
set.seed(647) # for reproducibility
porc_parental <- 0.4
idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
YTrn <- Y; YTrn[idPar,] <- NA

# ModelFit using BGLR
ETA <- list(list(X = X, model = "BRR"), list(K = K, model = "RKHS"))
model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,
                    resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = FALSE),
                    nIter = 30000, burnIn = 5000)

R <- as.matrix(model$resCov$R)               # residual cov matrix
B0 <- as.numeric(model$mu)                   # Overall mean
yHat <- model$ETAHat[model$missing_records,] # Puntual BVs

# Eval loss functions
out <- ApproxMPS(B0 = B0, yHat = yHat, R = R, method = "kl")

# Best 15% = 36 lines
selected_lines <- which(order(out$loss) %in% 1:36)
print(selected_lines)
out$yHat[selected_lines,]
```

**Example 3: Multi-trait selection with positive and negative direction of genetic progress**

```r
# Clean and setting local environment
rm(list = ls())
setwd("your path")
library(MPS); library(BGLR)

# Loading dataset
data(AdvEYT)
Xscaled <- scale(X, center = TRUE, scale = TRUE)
G <- tcrossprod(Xscaled) /ncol(X)
Y2 <- Y

direction <- c(1,1,1,1,1,1,1,1,-1,-1,-1)
Y2 <- sweep(Y2, 2, direction, '*') 

ETA <- list(list(K = G, model = "RKHS"))
model <- Multitrait(y = Y2, ETA = ETA, intercept = TRUE,
                    resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = FALSE),
                    nIter = 100000, burnIn = 20000,
                    saveAt = "Fit_")

# Evaluation of loss function
out <- ApproxMPS(B0 = as.numeric(model$mu), 
                 yHat = model$ETAHat , 
                 R = as.matrix(model$resCov$R), 
                 method = "kl",
                 direction = direction)

# Best 15% = 30 lines
selected_lines <- which(order(out$loss) %in% 1:30)
print(selected_lines)
sweep(out$yHat[selected_lines,], 2, direction, '*') 
```
