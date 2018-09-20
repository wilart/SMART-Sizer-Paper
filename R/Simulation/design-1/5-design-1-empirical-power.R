#Simulation Design 1: Compute the empirical power
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/Simulation/compute-best-set-empirical-power.R")
#run scripts 1-design-1-simulation.R, 2-design-1-empirical-covariance.R, 3-design-1-c-quantiles.R, and 4-design-1-compute-theta-AIPW.R first
n_grid <- seq(100, 500, 50)
in.best.set.indicator <- vector('list', length = length(n_grid))
load("~/../Dropbox/SMART-R/Rda/design.1.c.by.n.delta.0.25.rda")
for (i in 1:length(n_grid)){
  for (j in 1:1000){
    in.best.set.indicator[[i]][[j]] <- ComputeBestSet(theta = design.1.theta.by.n.delta.0.25[[i]][[j]], c_alpha = design.1.c.by.n.delta.0.25[[i]][[j]], V = design.1.cov.by.n.delta.0.25[[i]][[j]], sample_size = n_grid[i])
  }
}

design.1.empirical.power.by.n.delta.0.25 <- vector("numeric", length = length(n_grid))

for (i in 1:length(n_grid)){
  design.1.empirical.power.by.n.delta.0.25[i] <- ComputeEmpiricalPower(Delta = 0.5, theta.true, in.best.set.indicator.one = in.best.set.indicator[[i]])
}


#save(design.1.empirical.power.by.n.delta.0.25, file = "~/../Dropbox/SMART-R/Rda/design.1.empirical.power.by.n.delta.0.25.rda")


