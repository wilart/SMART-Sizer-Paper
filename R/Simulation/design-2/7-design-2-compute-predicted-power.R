#Simulation Design 2: Compute the predicted power over a grid of n-values

source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")

n_grid <- seq(100,500,50)

set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)

Delta <- c(2.752, 0.75, 1, 0, 0.75)

#Identity
## Covariance = Identity
V <- diag(5)
c_alpha <- computeC(V, mc_list_SVD = mc_list, alpha = 0.05) #indep. of n

design.2.predicted.power.identity <- rep(0, length(n_grid))

computeSampleSize(V, Delta = c(2.752, 0.75, 1, 0, 0.75), min_Delta = 0.7, desired_power = 0.8, c_alpha, Z = mc_list)


for (i in 1:length(n_grid)) {
  design.2.predicted.power.identity[i] <- computePower(V, Delta, min_Delta = 0.7, sample_size = n_grid[i], c_alpha, Z = mc_list)
}

###True
load('~/../Dropbox/SMART-R/Rda/design.2.cov.true.list.rda')
design.2.cov.true <- apply(simplify2array(design.2.cov.true.list),1:2, mean)

V <- design.2.cov.true
mc_list_SVD <- ComputeSVD(V, mc_list)
c_alpha <- computeC(V, mc_list_SVD, alpha = 0.05)
design.2.predicted.power.true <- rep(0, length(n_grid))

computeSampleSize(V, Delta = c(2.752, 0.75, 1, 0, 0.75), min_Delta = 0.7, desired_power = 0.8, c_alpha, Z = mc_list_SVD)


for (i in 1:length(n_grid)) {
  design.2.predicted.power.true[i] <- computePower(V, Delta, min_Delta = 0.7, sample_size = n_grid[i], c_alpha, Z = mc_list_SVD)
}


### No correlation conservative
V <-diag(diag(design.2.cov.true),5)

mc_list_SVD <- ComputeSVD(V, mc_list)
c_alpha <- computeC(V, mc_list_SVD, alpha = 0.05)


computeSampleSize(V, Delta = c(2.752, 0.75, 1, 0, 0.75), min_Delta = 0.7, desired_power = 0.8, c_alpha, Z = mc_list_SVD)


design.2.predicted.power.conservative <- rep(0, length(n_grid))

for (i in 1:length(n_grid)) {
  design.2.predicted.power.conservative[i] <- computePower(V, Delta, min_Delta = 0.7, sample_size = n_grid[i], c_alpha, Z = mc_list_SVD)
}

