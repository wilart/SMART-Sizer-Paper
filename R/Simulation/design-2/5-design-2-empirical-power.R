#Simulation Design 2: compute empirical power
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/Simulation/compute-best-set-empirical-power.R")


#run scripts 1-design-2-simulation.R, 2-design-2-empirical-covariance.R, 3-design-2-c-quantiles.R, and 4-design-2-compute-theta-AIPW.R first
load('~/../Dropbox/SMART-R/Rda/design.2.c.by.n.delta.2.rda')

n_grid <- seq(100, 500, 50)


in.best.set.indicator <- vector('list', length = length(n_grid))


for (i in 1:length(n_grid)){
  for (j in 1:1000){
    in.best.set.indicator[[i]][[j]] <- ComputeBestSet(theta = design.2.theta.by.n.delta.2[[i]][[j]], c_alpha = design.2.c.by.n.delta.2[[i]][[j]], V = design.2.cov.by.n.delta.2[[i]][[j]], sample_size = n_grid[i])
  }
}

design.2.empirical.power.by.n.delta.2 <- vector("numeric", length = length(n_grid))

for (i in 1:length(n_grid)){
  design.2.empirical.power.by.n.delta.2[i] <- ComputeEmpiricalPower(Delta = 0.7, design.2.theta.true, in.best.set.indicator.one = in.best.set.indicator[[i]])
}

#save(design.2.empirical.power.by.n.delta.2, file = '~/../Dropbox/SMART-R/Rda/design.2.empirical.power.by.n.delta.2.rda')
