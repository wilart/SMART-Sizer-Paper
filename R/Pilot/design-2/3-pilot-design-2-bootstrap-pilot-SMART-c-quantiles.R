#Compute the quantiles c
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")

set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.2.pilot.c.unknown.var <- foreach(k=1:1000) %dopar% {
  computeC(V = design.2.pilot.cov.bootstrap.list[[k]],
           mc_list_SVD = ComputeSVD(V = design.2.pilot.cov.bootstrap.list[[k]],mc_list), 
           alpha = 0.05)
}

stopCluster(cl)

#save(design.2.pilot.c.unknown.var, file = "~/../Dropbox/SMART-R/Rda/design.2.pilot.c.unknown.var.rda")


no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.2.pilot.c.known.var <- foreach(k=1:1000) %dopar% {computeC(V = design.2.pilot.cov.bootstrap.known.var[[k]],
                                                                       mc_list_SVD = ComputeSVD(V = design.2.pilot.cov.bootstrap.known.var[[k]], mc_list), 
                                                                       alpha = 0.05)
}

stopCluster(cl)
#save(design.2.pilot.c.known.var, file = "~/../Dropbox/SMART-R/Rda/design.2.pilot.c.known.var.rda")
