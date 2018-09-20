#Compute the quantiles c
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")

set.seed(1248789)
mc_list <- simulateNormal(n_EDTR = 4, n_sim = 1000)

library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.1.pilot.c.unknown.var <- foreach(k=1:1000) %dopar% {
  computeC(V = design.1.pilot.cov.bootstrap.list[[k]],
          mc_list_SVD = ComputeSVD(V = design.1.pilot.cov.bootstrap.list[[k]], mc_list), 
          alpha = 0.05)
}

stopCluster(cl)

#save(design.1.pilot.c.unknown.var, file = "~/../Dropbox/SMART-R/Rda/design.1.pilot.c.unknown.var.rda")
load("~/../Dropbox/SMART-R/Rda/design.1.pilot.c.unknown.var.rda")

no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)
design.1.pilot.c.known.var.list <- foreach(k=1:1000) %dopar% {computeC(V = design.1.pilot.cov.bootstrap.known.var[[k]],
                                      mc_list_SVD = ComputeSVD(V = design.1.pilot.cov.bootstrap.known.var[[k]], mc_list), 
                                      alpha = 0.05)
  }

stopCluster(cl)






#save(design.1.pilot.c.known.var.list, file = "~/../Dropbox/SMART-R/Rda/design.1.pilot.c.known.var.list.rda")
load("~/../Dropbox/SMART-R/Rda/design.1.pilot.c.known.var.list.rda")
#save(design.1.pilot.c, file = '~/../Dropbox/SMART-R/Rda/design.1.pilot.c.list.7.16.rda')
#save(design.1.c.by.n.delta.0.25, file = "C:/Users/Willi/Dropbox/SMART-R/Power-Analysis-in-a-SMART-design/rda/design.1.c.by.n.delta.0.25.rda")

#load("C:/Users/Willi/Dropbox/SMART-R/Power-Analysis-in-a-SMART-design/rda/design.1.c.by.n.delta.0.25.rda")
