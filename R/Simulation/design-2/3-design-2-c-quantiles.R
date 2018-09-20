#Simulation Design 2: Compute quantiles c
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")

set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)

design.2.c.by.n.delta.2 <- vector("list", length = length(n_grid))

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

design.2.c.by.n.delta.2 <- foreach(i = 1:length(n_grid)) %:%
  foreach(j = 1:n_sim) %dopar% {
    computeC(V = design.2.cov.by.n.delta.2[[i]][[j]],
             mc_list_SVD = ComputeSVD(V = design.2.cov.by.n.delta.2[[i]][[j]], mc_list), 
             alpha = 0.05)
  }

stopCluster(cl)


#save(design.2.c.by.n.delta.2, file = "~/../Dropbox/SMART-R/Rda/design.2.c.by.n.delta.2.rda")

