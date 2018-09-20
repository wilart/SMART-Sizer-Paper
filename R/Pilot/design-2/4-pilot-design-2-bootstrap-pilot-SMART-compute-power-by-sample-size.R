#Compute the bootstrap power confidence interval
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")

Delta.DR <- c(2.752, 0.75, 1, 0, 0.75)
n_grid <- seq(100, 600, 50)

#Unknown variances
design.2.pilot.power.list <- vector("list", length = length(n_grid))
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.2.pilot.power.list <- foreach(j = 1:length(n_grid)) %:%
  foreach(i = 1:1000) %dopar% {
    computePower(V = design.2.pilot.cov.bootstrap.list[[i]], Delta = Delta.DR, min_Delta = 0.7, sample_size = n_grid[j], c_alpha = design.2.pilot.c.unknown.var[[i]], Z = ComputeSVD(V = design.2.pilot.cov.bootstrap.list[[i]], mc_list))
  }

stopCluster(cl)

#save(design.2.pilot.power.list, file = "~/../Dropbox/SMART-R/Rda/design.2.pilot.power.list.rda")



no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.2.pilot.power.known.var.list <- foreach(j = 1:length(n_grid)) %:%
  foreach(i = 1:1000) %dopar% {
    computePower(V = design.2.pilot.cov.bootstrap.known.var[[i]], Delta = Delta.DR, min_Delta = 0.7, sample_size = n_grid[j], c_alpha = design.2.pilot.c.known.var[[i]], Z = ComputeSVD(design.2.pilot.cov.bootstrap.known.var[[i]], mc_list))
  }

stopCluster(cl)
#save(design.2.pilot.power.known.var.list, file = "~/../Dropbox/SMART-R/Rda/design.2.pilot.power.known.var.list.rda")
load("~/../Dropbox/SMART-R/Rda/design.2.pilot.power.known.var.list.rda")
load("~/../Dropbox/SMART-R/Rda/design.2.pilot.power.list.rda")

load('~/../Dropbox/SMART-R/Rda/design.2.cov.true.list.rda')
design.2.cov.true <- apply(simplify2array(design.2.cov.true.list),1:2, mean)


set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)

c_alpha_true = computeC(design.2.cov.true, ComputeSVD(design.2.cov.true, mc_list), alpha = 0.05)
c_alpha_conservative = computeC(diag(diag(design.2.cov.true),5), ComputeSVD(diag(diag(design.2.cov.true),5), mc_list), alpha = 0.05)
c_sample.known.var <- computeC(V= design.2.pilot.sample.known.var.cov, mc_list_SVD = ComputeSVD(V = design.2.pilot.sample.known.var.cov, mc_list), alpha = 0.05)
c_sample.unknown.var <- computeC(V= design.2.pilot.sample.cov, mc_list_SVD = ComputeSVD(V = design.2.pilot.sample.cov, mc_list), alpha = 0.05)


pilot.power.summary <- cbind(n_grid, "True Power" = sapply(n_grid, function(x) computePower(V = design.2.cov.true, Delta.DR , min_Delta = 0.7, c_alpha = c_alpha_true, sample_size = x, Z = ComputeSVD(design.2.cov.true, mc_list))),
                             "Pilot Unknown Sigma" = sapply(n_grid, function(x) computePower(V = design.2.pilot.sample.cov, Delta.DR, min_Delta = 0.7, sample_size = x, c_alpha = c_sample.unknown.var, Z = ComputeSVD(design.2.pilot.sample.cov, mc_list))),
                             "Pilot Known Sigma" = sapply(n_grid, function(x) computePower(V = design.2.pilot.sample.known.var.cov, Delta.DR, min_Delta = 0.7, sample_size = x, c_alpha = c_sample.known.var, Z = ComputeSVD(design.2.pilot.sample.known.var.cov, mc_list))),
                             "Conservative Power" = sapply(n_grid, function(x) computePower(V = diag(diag(design.2.cov.true),5), Delta.DR , min_Delta = 0.7, c_alpha = c_alpha_conservative, sample_size = x, Z = ComputeSVD(diag(diag(design.2.cov.true),5), mc_list))))

pilot.known.var.CI <- cbind(n_grid,do.call(rbind, lapply(design.2.pilot.power.known.var.list, function(x) quantile(unlist(x), probs = c(0.025,0.975)))))
pilot.unknown.var.CI <- cbind(n_grid,do.call(rbind, lapply(design.2.pilot.power.list, function(x) quantile(unlist(x), probs = c(0.025,0.975)))))
