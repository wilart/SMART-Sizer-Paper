#Compute the bootstrapped power confidence interval
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")

Delta.DR <- c(0, 0.502, 0.103, 0.605)
n_grid <- seq(100, 600, 50)

#Unknown variances
design.1.pilot.unknown.var.power.list <- vector("list", length = length(n_grid))
library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.1.pilot.unknown.var.power.list <- foreach(j = 1:length(n_grid)) %:%
  foreach(i = 1:1000) %dopar% {
    computePower(V = design.1.pilot.cov.bootstrap.list[[i]], Delta = Delta.DR, min_Delta = 0.5, sample_size = n_grid[j], c_alpha = design.1.pilot.c.unknown.var[[i]], Z = ComputeSVD(design.1.pilot.cov.bootstrap.list[[i]],mc_list))
}

stopCluster(cl)

#save(design.1.pilot.unknown.var.power.list, file = "~/../Dropbox/SMART-R/Rda/design.1.pilot.unknown.var.power.list.rda")


no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.1.pilot.power.known.var.list <- foreach(j = 1:length(n_grid)) %:%
  foreach(i = 1:1000) %dopar% {
    computePower(V = design.1.pilot.cov.bootstrap.known.var[[i]], Delta = Delta.DR, min_Delta = 0.5, sample_size = n_grid[j], c_alpha = design.1.pilot.c.known.var.list[[i]], Z = ComputeSVD(design.1.pilot.cov.bootstrap.known.var[[i]], mc_list))
  }

stopCluster(cl)
#save(design.1.pilot.power.known.var.list, file = "~/../Dropbox/SMART-R/Rda/design.1.pilot.power.known.var.list.rda")
load("~/../Dropbox/SMART-R/Rda/design.1.pilot.power.known.var.list.rda")
load("~/../Dropbox/SMART-R/Rda/design.1.pilot.unknown.var.power.list.rda")

library(Matrix)
load("~/../Dropbox/SMART-R/Rda/design.1.cov.true.list.rda")

design.1.cov.true <- nearPD(apply(simplify2array(design.1.cov.true.list),1:2, mean))$mat

set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 4, n_sim = 1000)

c_alpha_true = computeC(design.1.cov.true, ComputeSVD(design.1.cov.true, mc_list), alpha = 0.05)
c_alpha_conservative = computeC(diag(diag(design.1.cov.true),4), ComputeSVD(diag(diag(design.1.cov.true),4), mc_list), alpha = 0.05)
c_sample.known.var <- computeC(V = design.1.pilot.sample.known.var.cov, mc_list_SVD = ComputeSVD(V = design.1.pilot.sample.known.var.cov, mc_list), alpha = 0.05)
c_sample.unknown.var <- computeC(V = design.1.pilot.sample.cov, mc_list_SVD = ComputeSVD(V = design.1.pilot.sample.cov, mc_list), alpha = 0.05)


pilot.power.summary <- cbind(n_grid, "True Power" = sapply(n_grid, function(x) computePower(V = design.1.cov.true, Delta.DR , min_Delta = 0.5, c_alpha = c_alpha_true, sample_size = x, Z = ComputeSVD(design.1.cov.true, mc_list))),
      "Pilot Unknown Sigma" = sapply(n_grid, function(x) computePower(V = design.1.pilot.sample.cov, Delta.DR, min_Delta = 0.5, sample_size = x, c_alpha = c_sample.unknown.var, Z = ComputeSVD(design.1.pilot.sample.cov, mc_list))),
      "Pilot Known Sigma" = sapply(n_grid, function(x) computePower(V = design.1.pilot.sample.known.var.cov, Delta.DR, min_Delta = 0.5, sample_size = x, c_alpha = c_sample.known.var, Z = ComputeSVD(design.1.pilot.sample.known.var.cov, mc_list))),
      "Conservative Power" = sapply(n_grid, function(x) computePower(V = diag(diag(design.1.cov.true),4), Delta.DR , min_Delta = 0.5, c_alpha = c_alpha_conservative, sample_size = x, Z = ComputeSVD(diag(diag(design.1.cov.true),4), mc_list))))

pilot.known.var.CI <- cbind(n_grid,do.call(rbind, lapply(design.1.pilot.power.known.var.list, function(x) quantile(unlist(x), probs = c(0.025,0.975)))))
pilot.unknown.var.CI <- cbind(n_grid,do.call(rbind, lapply(design.1.pilot.unknown.var.power.list, function(x) quantile(unlist(x), probs = c(0.025,0.975)))))
