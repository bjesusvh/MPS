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

# Best 15% = 36 lines
selected_lines <- which(order(out$loss) %in% 1:36)
print(selected_lines)
out$yHat[selected_lines,]
```

**Example 2: Multi-trait Genomic + pedigree selection**

```r
# Cleaning and setting local environment
rm(list = ls())
setwd("~/Desktop/R Package/Ejemplos Finales/example2")
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
