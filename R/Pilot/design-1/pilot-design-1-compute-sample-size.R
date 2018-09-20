#Compute sample size for the algorithm given in Section 7 to compute bootstrap maximum sample size
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")

set.seed(2374)
mc_list <- simulateNormal(n_EDTR = 4, n_sim = 1000)
Delta.DR <- c(0, 0.502, 0.103, 0.605)



library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

sample_size.bootstrap.list <- foreach(i = 1:1000) %dopar% {
  computeSampleSize(V= design.1.pilot.cov.bootstrap.known.var[[i]],Delta = Delta.DR, min_Delta = 0.5, desired_power = 0.8, c_alpha = computeC(V = design.1.pilot.cov.bootstrap.known.var[[i]], mc_list_SVD = ComputeSVD(V = design.1.pilot.cov.bootstrap.known.var[[i]], mc_list)),Z = ComputeSVD(V= design.1.pilot.cov.bootstrap.known.var[[i]], mc_list))
}


stopCluster(cl)
sample_size.bootstrap.list
max(unlist(sample_size.bootstrap.list))
design.1.pilot.cov.bootstrap.known.var
load("~/../Dropbox/SMART-R/Rda/sample_size.bootstrap.list.rda")
#save(sample_size.bootstrap.list, file = "~/../Dropbox/SMART-R/Rda/sample_size.bootstrap.list.rda")
