# MPS
An R Package to Assist Multi-trait Parental Selection via Genomic, Pedigree, or Phenotypic Selection.


## Examples

**Example using Genomic Selection**

```{r}
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
set.seed(384)
porc_parental <- 0.4
idPar <- sample(1:n, ceiling(porc_parental*n), replace = FALSE)
XTrn <- X[-idPar,]
YTrn <- Y[-idPar,]

# ModelFit using BGLR
ETA <- list(list(X = XTrn, model = "BRR", saveEffects = TRUE))
model <- Multitrait(y = YTrn,
  ETA = ETA,
  intercept = TRUE,
  resCov = list(df0 = 5, S0 = NULL,type = "UN", saveEffects = TRUE),
  saveAt = paste0("./chains/", "BRR_"),
  nIter = 100000,
  burnIn = 30000)

# Reading Posterior MCMC
mu <- as.matrix(read.table(file = "./chains/BRR_mu.dat", header = FALSE))
B <- readBinMatMultitrait('./chains/BRR_ETA_1_beta.bin')
R <- as.matrix(read.table(file = "./chains/BRR_R.dat", header = FALSE))
XPar <- X[idPar, ]


# Evaluating Loss Function
out <- FastMPS(Xcand = XPar,
          B0 = mu,
          B = B,
          R = R,
          p = 0.1,
          method = "energy")

# Plotting results
colnames(out$yHat) <- colnames(Y)
pairs(out$yHat,
  col = ifelse(out$selected, "#2c7fb8", "darkgray"),
  pch = ifelse(out$selected, 19, 1))
```


**Example using Genomic + Pedigree**

**Example using Phenotypic selection**

**Example for Multi-environment best parent selection**
