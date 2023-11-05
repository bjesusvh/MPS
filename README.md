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

In this example, we use the wheat database that is included with the `R` package `BGLR`. The database consists of genotypic and phenotypic information of 599 wheat lines in four traits. The genotypic information comprises of 1279 SNP-type molecular markers. In the following code, the `R` memory is cleaned, and the working directory is fixed. After, the two R packages used are loaded into `R` (`BGLR` for regression model, `MPS` for evaluation of posterior expected loss). Phenotypic records are saved on `Y` object and `X` contains the standardized genotypic data (SNPs molecular markers).

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
```
After that, the data is partitioned, with 60% being used to train the model for prediction purposes (simulating the scenario of having related historical information), while the remaining 40% of the data simulates the set of candidate individuals for selection. set.seed(647) is used only for reproducibility of results presented in this article. The ID of candidate/parental lines is sampled without replacement and then sorted and saved in idPar R vector. The last three lines of code perform a sub setting of phenotypic and genotypic data in training and parental set.

```r
XTrn <- X[-idPar,]        # Genotypic data in Training set
XPar <- X[idPar, ]        # Genotypic data in Parental set
YTrn <- Y[-idPar,]        # Phenotypic data in Training set
```
Next, the linear predictor is specified in ETA object imposing Gaussian priors (BRR) on regression coefficients to emulate Ridge Regression. The multi-trait model is fitted using the Multi-trait function of BGLR R package passing as arguments the phenotypic records and the linear predictor. Note that saveEffects = TRUE in the linear predictor specification (ETA) and in resCov argument is for saving posterior MCMC samples of regression coefficients and residual variance-covariance, respectively.

```r
# Specification of the linear predictor for BGLR
ETA <- list(list(X = X, model = "BRR", saveEffects = TRUE))

# Run the model using Multitrait function of BGLR R Package
model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,   
                    resCov = list(type = "UN", saveEffects = TRUE), 
                    nIter = 100000, burnIn = 30000)
```
The following lines of code illustrate how to read the posterior samples from the Markov chain Monte Carlo (MCMC) of the model’s parameters. It is necessary to verify that the MCMC samples have converged to their stationary distributions (CODA package of Plummer (2006) can assist). After the aforementioned step, we pass these MCMC samples as input to the FastMPS function. Next, we will explain each input in detail.
In this example, `Xcand = XPar` correspond to the matrix of predictors (scaled incidence matrix based on SNPs) for the candidate lines of selection. The dimension of `XPar` is of `n=240` rows (lines) and `k=1279` columns (predictors). `B0` is a matrix of dimension `M=10000 × t=4` corresponding to the intercept/overall mean term in the linear model, where `M` is the number of MCMC samples and `t` is the number of traits. The next argument is `B` which is an array containing regression coefficients of dimension `M=10000 × k=1279 × t=4` of MCMC samples. `R` input represents a matrix of dimension `M×(t×(t+1)/2)` of Markov Chain Monte Carlo samples of the variance-covariance components in the residual covariance matrix. In this example `R` matrix has dimensions `10000 × 10`. In the last argument of FastMPS function is fixed as `"kl"` for Kullback-Leibler Loss. The last line of code obtains the average similarity between each line with the remaining candidates.


```r
# Reading Posterior MCMC of model parameters to compute expected loss
B0 <- as.matrix(read.table(file = "mu.dat",  # Overall mean
                           header = FALSE))
B <- readBinMatMultitrait('ETA_1_beta.bin')  # Regression coefficients
R <- as.matrix(read.table(file = "R.dat",    # Residual covariance matrix
                          header = FALSE))

# Evaluate loss function
out <- FastMPS(Xcand = XPar, B0 = B0, B = B, R = R, method = "kl")
out
```
After running the FastMPS function, the output list contains the posterior expected loss (loss), the ranking of each candidate line for selection, and a data.frame (yHat) containing the predictions (BV) for all traits.

```r
> str(out)
List of 4
 $ method : chr "kl"
 $ loss   : num [1:240] 0.482 0.466 0.222 0.182 0.49 ...
 $ ranking: num [1:240] 214 208 71 41 217 236 201 59 126 158 ...
 $ yHat   :'data.frame':    240 obs. of  4 variables:
  ..$ X1: num [1:240] -0.833 -0.847 0.497 0.511 -0.125 ...
  ..$ X2: num [1:240] 0.186 0.139 -0.394 0.211 -0.348 ...
  ..$ X3: num [1:240] -0.366 -0.4283 -0.2819 0.1197 0.0116 ...
  ..$ X4: num [1:240] -0.276 -0.2762 -0.3616 -0.0803 0.0753 ...
 - attr(*, "class")= chr "MPS"
```

**Example 2: Multi-trait Genomic + pedigree selection**

For this example, the same database as in the previous exercise (example 1) is employed, but with the incorporation of pedigree and molecular markers as predictors. The first part of the code is very similar to that of example 1, so we omit its explanation. The only difference here is that in K, the pedigree information is stored.

```r
rm(list = ls())                  # Cleaning memory
setwd("path in your computer")   # Working directory  
library(MPS); library(BGLR)      # Needed R Packages
data(wheat)                      # Loading dataset

Y <- as.matrix(wheat.Y)          # Phenotypic records  
X <- scale(wheat.X, center = TRUE, scale = TRUE) # Genotypic data
K <- wheat.A  
```
In total 60% of the data is allocated for training the model, while the remaining 40% simulates the role of candidate lines for selection. 

```r
set.seed(647)             # Seed for reproducibility
n <- nrow(Y)              # Number of records
porc_parental <- 0.4      # Percentage in candidates
idPar <- sort(sample(1:n, # Random observation for candidate set
                     ceiling(porc_parental*n), 
                     replace = FALSE))
YTrn <- Y; YTrn[idPar,] <- NA
```
The linear predictor now incorporates the matrix of molecular markers, and associated regression coefficients are subject to Ridge penalization. Meanwhile, the relationship matrix (pedigree) is used to incorporate random effects which is modeled through nonparametric genomic regressions based on reproducing kernel Hilbert spaces (RKHS) methods, details of RKHS approach can be found in Perez-Rodriguez and de los Campos (2022).


```r
# ModelFit using BGLR
ETA <- list(list(X = X, model = "BRR"), list(K = K, model = "RKHS"))
model <- Multitrait(y = YTrn, ETA = ETA, intercept = TRUE,
        resCov = list(type = "UN"), nIter = 100000, burnIn = 30000)
```
Finally, outputs from Multi-trait R function are inputs of the MPS R Package. In this example, the ApproxMPS function is used, which is a version that approximates the loss of each candidate line for selection when the MCMC samples are not available from BGLR. This approximation should be taken with caution, as it assumes that the predictive distribution of each candidate line is multivariate normal.
The ApproxMPS function receives several arguments. Firstly, B0 is a vector representing the overall mean of each trait, with a length equal to the number of traits (t=4). Then, yHat should be a matrix containing punctual posterior mean of each candidate, with dimensions n=240 x t=4, where n is the number of candidates and t is the number of traits. The next argument is R, which is a square matrix representing the residual covariance matrix with dimensions t=4 x t=4.
The target argument is an optional vector of length equal to the number of traits (t=4) that reflects the breeder’s expectation. If not provided, the algorithm determines it internally. Lastly, the method argument specifies the loss function to be used and must be one of "kl" for Kullback-Leibler, "energy" for Energy Score, or "malf" for Multivariate Asymmetric Loss.

```r
B0 <- as.numeric(model$mu)                    # Overall mean
yHat <- model$ETAHat[model$missing_records, ] # Puntual BVs
R <- as.matrix(model$resCov$R)                # Residual Covariance matrix
out <- ApproxMPS(B0 = B0, yHat = yHat,        # Evaluate loss functions
     R = R, method = "kl") 

# Best 15% = 36 lines
selected_lines <- which(order(out$loss) %in% 1:36)
print(selected_lines)
```

**Example 3: Multi-trait selection with positive and negative direction of genetic progress**

In this example we employed a database named AdvEYT which is included in the MPS R package. The data comprised phenotypic records for 190 lines for eleven traits: GY_B5IR_F5IR_BEHT, GY_BLHT, GY_B2IR, and GY_TPE belonging to the grain yield category; Zinc, TKW, GRNPRO and ALVW fall under the quality category; and YR_LUD, YR_NJ and SR_NJ correspond to the diseases category. Furthermore, the dataset includes information on 8,545 SNP-type molecular markers.
The objective of this analysis is to identify the top 15 percent of lines (approximately 28 lines) based on the criterion of minimum expected loss. For the Grain Yield and Quality traits, the desired direction of genetic progress is positive, whereas for the Diseases category traits, progress in the negative direction is sought.
The first code segment, clears the R memory, sets the desired working directory, loads the necessary libraries and the database to be used. For computational convenience we used the genomic relationship matrix (G=(ZZ^')/p) to fit a multi-trait GBLUP model, where Z is the standardized SNP matrix. The last line of code creates a copy of the phenotypic records matrix to reverse direction of observed phenotypic values for Disease category traits.

```r
rm(list = ls());                                 # Clean R memory
setwd("your-path")                               # Set the working directory    
library(MPS); library(BGLR)                      # Load R Packages
data(AdvEYT)                                     # Load data
Z <- scale(X, center = TRUE, scale = TRUE)       # Standardized SNPs
G <- tcrossprod(Z) /ncol(Z)                      # G relationship matrix
Y2 <- Y                                          # copy of the response                  
```
The next code begins by defining a vector that indicates the desired direction of improvement: 1 indicates that we aim to increase the genetic value in that trait, while -1 indicates a desire to decrease the genetic value. In this case, the phenotypic records of the traits YR_LUD, YR_NJ, and SR_NJ are located in the last three columns of Y, therefore direction object has eight ones and finally three minus one. Subsequently, in Y2, we reverse the direction of the response variable. We define the linear predictor (ETA) as a function of the G matrix by incorporating random effects, which are modeled through nonparametric genomic regressions based on reproducing kernel Hilbert spaces (RKHS) methods. This approach allows BGLR to achieve computationally efficient model fitting, however, it does not yield Markov Chain Monte Carlo chains. Therefore, we will use the ApproxMPS function from the MPR package to approximate expected losses. Note that Y2 is given to the regression model using the Multitrait function from BGLR, not Y.

```r
direction <- c(1,1,1,1,1,1,1,1,-1,-1,-1)
Y2 <- sweep(Y2, 2, direction, '*') 
ETA <- list(list(K = G, model = "RKHS"))
model <- Multitrait(y = Y2, ETA = ETA, intercept = TRUE,
          resCov = list(type = "UN"),
          nIter = 100000, burnIn = 20000)
```
Finally, the ApproxMPS function is called, and arguments like example 2 are provided, which are the outputs of BGLR. 

```r
out <- ApproxMPS(B0 = as.numeric(model$mu), yHat = model$ETAHat , 
                 R = as.matrix(model$resCov$R), method = "kl",
                 direction = direction)

# Best 15% = 30 lines
selected_lines <- which(order(out$loss) %in% 1:30)
print(selected_lines)
sweep(out$yHat[selected_lines,], 2, direction, '*') 
```

**Example 5. Candidates’ selection with missing phenotypic records for some traits**

At times, partial phenotypic information is available for certain traits in some candidate lines for selection. By combining prediction using BGLR and expected loss calculation using MPS, it is possible to identify which lines the breeder should select. To illustrate this concept, the database accompanying the MPS package is used. Let’s assume we are interested in identifying the top 20% (40 lines) considering the traits GY_BLHT, TKW, GRNPRO, and SR_NJ. However, out of the 190 lines in the dataset, 70 of them lack phenotypic records for the trait GY_BLHT, 50 lines have no information on TKW, 60 lines have missing data for GRNPRO, and 20 lines have missing information for SR_NJ.
The following set of code clears the R memory, sets the working directory, loads the necessary packages, generates a seed for reproducibility in this example, and subsets the relevant traits.


```r
rm(list = ls())                                        # Clean Work space
setwd("~/Desktop/R Package/Ejemplos Finales/example5") # Working directory
library(MPS); library(BGLR)                            # Needed libraries
set.seed(752); data(AdvEYT)                       # Reproducibility and Load dataset
Y <- Y[,c('GY_BLHT', 'TKW', 'GRNPRO', 'SR_NJ')]   # Traits of interest
n <- nrow(Y)                                      # Number of candidates
```
To simulate the scenario of missing data, randomly missing data (NA) are introduced. As improvement direction is decreasing for the third trait, it is multiplied by -1.

```r
# Allocating randomly missing phenotypic data in each trait
Y[sort(sample(1:n, size = 70, replace = FALSE)), 1] <- NA
Y[sort(sample(1:n, size = 50, replace = FALSE)), 2] <- NA
Y[sort(sample(1:n, size = 60, replace = FALSE)), 3] <- NA
Y[sort(sample(1:n, size = 20, replace = FALSE)), 4] <- NA
Y2 <- Y                   # Copy of phenotypic records
Y2[,4] <- Y2[,4] * -1     # Change sign in SR_NJ trait
```

The genotypic information is standardized, and the linear predictor (ETA) is defined by selecting model = SpikeSlab to instruct BGLR to impose prior distributions on regression coefficients, allowing for automatic variable selection (for details, see Pérez-Rodrı́guez & de los Campos, 2022). The argument saveEffects = TRUE instructs to save the MCMC chains. Finally, the model is fitted by specifying the number of iterations (nIter) and the burn-in (burnIn). The argument resCov = list(type = "UN", saveEffects = TRUE) indicates to BGLR that an unstructured covariance matrix is assumed for the residual covariance matrix of the linear model, and it is also requested to save the MCMC chains for this matrix.

```r
# Center and scaling genomic info
X <- scale(X, center = TRUE, scale = TRUE)

# Specification of the linear predictor for BGLR
ETA <- list(list(X = X, model = "SpikeSlab", saveEffects = TRUE))

# Run the model
model <- Multitrait(y = Y2, ETA = ETA, intercept = TRUE,
 resCov = list(type = "UN", saveEffects = TRUE), nIter = 100000, burnIn = 30000)
```
To obtain the expected loss, it is necessary to read the MCMC chains from the model fit, which serves as input to the FastMPS function. The direction argument specifies that for the first three traits, the genetic value need to be increased, while for the fourth trait, the progress is in the opposite direction. Finally, the average values of the predicted BVs for the selected and non-selected lines are presented.

```r
# Reading Posterior MCMC of model parameters & direction
B0 <- as.matrix(read.table(file = "mu.dat", header = FALSE))# Overall mean                
B <- readBinMatMultitrait('ETA_1_beta.bin')                 # Regression coefficients
R <- as.matrix(read.table(file = "R.dat", header = FALSE))  # Residual cov. matrix
direction <- c(1,1,1,-1)                                   # Direction of improvement
out <- FastMPS(Xcand = X, B0 = B0, B = B, R = R,       # Evaluation of Loss Function
               method = "kl", direction = direction)

idToSelect <- out$ranking %in% 1:40            # ID of selected
yHat <- out$yHat                               # Predicted BVs
yHat[,4] <- -1*yHat[,4]                        # Reverse order of BVs for trait 4

colMeans(yHat[idToSelect,])                    # Average BVs of selected
       X1        X2        X3        X4 
101.85855  53.29620  13.35975  22.59924 

colMeans(yHat[!idToSelect,])                   # Average BVs of non-selected
       X1        X2        X3        X4 
100.95011  52.45389  13.10432  26.27064
```
The IDs of the selected individuals are presented below. The number in parentheses indicates the count of missing trait records. For instance, '1(3)' denotes that the top line (with the lowest expected loss) is the first candidate in the dataset, originally missing data in three traits. The same explanation applies to the rest of the selected candidates. This example illustrates how breeders can effectively choose individuals even in the presence of missing phenotypic records for some traits.

```r
NAsBySelLines <- apply(Y2[which(idToSelect),], 1, function(x) sum(is.na(x)))
IDSel <- which(idToSelect)
Print(paste0(IDSel, '(', NAsBySelLines, ')'))

[1] "1(3)"   "4(2)"   "14(1)"  "21(4)"  "29(1)"  "34(2)"  "40(2)"  "42(3)" 
[9] "49(1)"  "50(3)"  "53(1)"  "54(2)"  "58(3)"  "63(3)"  "65(2)"  "66(1)" 
[17] "70(3)"  "73(1)"  "79(2)"  "92(1)"  "95(0)"  "99(1)"  "101(3)" "118(2)"
[25] "119(1)" "129(0)" "131(3)" "137(1)" "142(3)" "143(2)" "144(2)" "147(1)"
[33] "148(3)" "149(1)" "160(1)" "162(2)" "171(2)" "173(2)" "174(3)" "190(2)"
```

**Example 6: Identifying superior candidates in multi-environment setting**

Sometimes, we have records for a specific trait assessed in different individuals and environments. MPS can also assist in identifying, from a set of candidates, which individuals perform best across all environments. In this example, we have taken data from 250 lines evaluated in three environments ('Bed5IR,' 'Drip,' and 'Bed2IR'). The original data, labeled 'DATASET4to5.Wheat.rdata,' can be accessed at https://hdl.handle.net/11529/10907, and its full description can be found in Montesinos-López et al. (2017). 

To borrow information across environments while allowing marker effects to change across environments, we fitted a Markers-Environment interaction model (MxE GBLUP) (López-Cruz et al., 2015). Detailed explanation how to fit MxE GBLUP model using BGLR can be found in https://github.com/MarcooLopez. To simulate scenarios of Training and Candidate/Parental sets, the data was randomly split into 20% (50 lines) for selection candidates. After fitting the model and approximating the expected loss for every candidate, we retained the top 40% of candidates (20 lines) with the minimum expected loss.

```r
rm(list = ls())                                        # Clean Work space
setwd("~/Desktop/R Package/Ejemplos Finales/example6") # Working directory
library(MPS); library(BGLR); library(tidyverse)        # Needed libraries
load('MultiEnvWheat.RData')                            # Load dataset                   
df <- Trigo.Data[,c(1,2,3)]                            # Keep only DTHD trait

PhenoEnv <- as.data.frame(pivot_wider(df, names_from = Env, values_from = DTHD))

# Seed for reproducibility
set.seed(647)                

# Geno and Pheno data
G <- as.matrix(G)

# Given records were adjusted by experimental design there are
# negative values, here we sum 100 units to put in positive values
Y <- PhenoEnv[,2:4] + 100
rownames(Y) <- PhenoEnv$Gid

# Simulating some NAs entries that correspond to candidates
idNA <- sample(x = 1:nrow(Y), size = nrow(Y)*0.2)
YNA <- Y
YNA[idNA, ] <- NA

# Design matrix for individuals. It connects individuals with environments
GID <- factor(rep(rownames(Y), ncol(Y)), levels=rownames(Y))
Zg <- model.matrix(~GID-1)

# Design matrix for environments. Used in the multi-environment R-Norm model
envID <- factor(rep(colnames(Y),each = nrow(Y)), levels = colnames(Y))
ZE <- model.matrix(~envID-1)   

#  Covariance structure for effects
ZgGZgt <- Zg %*% G %*% t(Zg)    # Genetic effect  
ZEZEt <- tcrossprod(ZE)         # Environmental effect

# Eigen decomposition (to speed computational time)
eigen_G <- eigen(G)
eigen_G0 <- eigen(ZgGZgt)

# Interaction terms (MxE model)
MxE_eigen <- vector("list",ncol(Y))
for(env in 1:ncol(Y)){ 
    tmp <- rep(0,ncol(Y)) ; tmp[env] <- 1; G1 <- kronecker(diag(tmp),G)
    MxE_eigen[[env]] <- eigen(G1)
}

n <- nrow(YNA)
nEnv <- ncol(YNA)
y <- as.numeric(unlist(YNA))

# Number of iterations and burn-in for Bayesian models
nIter <- 30000; burnIn <- 2000

# MxE interaction model
ETA <- list(list(~envID-1,model="FIXED"))
ETA[[2]] <- list(V=eigen_G0$vectors,d=eigen_G0$values,model='RKHS')

# Adding interaction terms
for(env in 1:nEnv){
    eigen_G1 <- MxE_eigen[[env]]
    ETA[[(env+2)]] <- list(V=eigen_G1$vectors,d=eigen_G1$values,model='RKHS')
}

# Model Fitting
fm <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)

# Retrieving point estimates
B0 = fm$mu + fm$ETA[[1]]$b
yHat = matrix(fm$yHat,ncol=nEnv)[idNA, ]
v1 = fm$varE + fm$ETA[[3]]$varU
v2 = fm$varE + fm$ETA[[4]]$varU
v3 = fm$varE + fm$ETA[[5]]$varU
R = diag(c(v1, v2, v3))

# Eval loss function
out <- ApproxMPS(B0 = B0, 
                 yHat = yHat, 
                 R = R, 
                 method = "kl")


# Only 40% of candidates is selected
idSel = out$ranking %in% 1:20

meansSel = colMeans(out$yHat[idSel,]) 
meansNoSel = colMeans(out$yHat[!idSel,]) 
dif1 = (meansSel - meansNoSel)
dif2 = (meansSel -  B0)
dif3 = (meansSel -  fm$mu )

dfToPlot = matrix(c(dif2, dif3), nrow = 2, byrow = TRUE)
colnames(dfToPlot) = colnames(Y)
rownames(dfToPlot) = c('With respect overall mean', 'With respect environment means')

barplot(dfToPlot, 
        main = "", 
        xlab = "Environments", 
        ylab = "Difference in Average BVs",
        beside = TRUE, 
        legend.text = rownames(dfToPlot),
        args.legend = list(bty = "n", x = "topleft", ncol = 1))
```
