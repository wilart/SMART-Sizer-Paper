#Simulation Design 1: Compute the predicted power over a grid of n-values

source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")

#Run design-1-covariance-true.R first
n_grid <- seq(100,500,50)
set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 4, n_sim = 1000)

Delta.DR <- c(0, 0.502, 0.103, 0.605)

## Covariance = Identity
V <- diag(4)


c_alpha <- computeC(V, mc_list_SVD = mc_list, alpha = 0.05) #indep. of n

computeSampleSize(V,Delta =c(0, 0.502, 0.103, 0.605), min_Delta = 0.5, desired_power = 0.8, c_alpha, Z = mc_list)


design.1.predicted.power.identity <- rep(0, length(n_grid))
for (i in 1:length(n_grid)) {
  design.1.predicted.power.identity[i] <- computePower(V, Delta = Delta.DR, min_Delta = 0.5, sample_size = n_grid[i], c_alpha, Z = mc_list)
}

#### Covariance = True: run design-1-covariance-true.R first
load("~/../Dropbox/SMART-R/Rda/design.1.cov.true.list.rda")
library(Matrix)
design.1.cov.true <- nearPD(apply(simplify2array(design.1.cov.true.list),1:2, mean))$mat
V <- design.1.cov.true
mc_list_SVD <- ComputeSVD(V, mc_list)
c_alpha <- computeC(V, mc_list_SVD)

computeSampleSize(V,Delta =c(0, 0.502, 0.103, 0.605), min_Delta = 0.5, desired_power = 0.8, c_alpha, Z = mc_list_SVD)


design.1.predicted.power.true <- rep(0, length(n_grid))

for (i in 1:length(n_grid)) {
  design.1.predicted.power.true[i] <- computePower(V, Delta = Delta.DR, min_Delta = 0.5, sample_size = n_grid[i], c_alpha, Z = mc_list_SVD)
}

####


#Conservative
V <- diag(diag(design.1.cov.true),4)
mc_list_SVD <- ComputeSVD(V, mc_list)
c_alpha <- computeC(V, mc_list_SVD, alpha = 0.05) #indep. of n
design.1.predicted.power.nocor <- rep(0, length(n_grid))

computeSampleSize(V,Delta =c(0, 0.502, 0.103, 0.605), min_Delta = 0.5, desired_power = 0.8, c_alpha, Z = mc_list_SVD)


for (i in 1:length(n_grid)) {
  design.1.predicted.power.nocor[i] <- computePower(V, Delta = Delta.DR, min_Delta = 0.5, sample_size = n_grid[i], c_alpha, Z = mc_list_SVD)
  
}
