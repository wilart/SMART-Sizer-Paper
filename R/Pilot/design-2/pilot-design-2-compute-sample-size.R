#Compute the bootstrap sample size and choose maximum sample size.
#Algorithm from Section 7 for choosing the sample size given a pilot SMART
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")

set.seed(2374)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)
Delta.DR <- c(2.752, 0.75, 1, 0, 0.75)



library(doParallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

sample_size.bootstrap.list.design.2 <- foreach(i = 1:1000) %dopar% {
  computeSampleSize(V= design.2.pilot.cov.bootstrap.known.var[[i]],Delta = Delta.DR, min_Delta = 0.7, desired_power = 0.8, c_alpha = computeC(V = design.2.pilot.cov.bootstrap.known.var[[i]], mc_list_SVD = ComputeSVD(V = design.2.pilot.cov.bootstrap.known.var[[i]], mc_list)),Z = ComputeSVD(V= design.2.pilot.cov.bootstrap.known.var[[i]], mc_list))
}


stopCluster(cl)

sample_size.bootstrap.list.design.2
max(unlist(sample_size.bootstrap.list.design.2))
design.2.pilot.cov.bootstrap.known.var
load("~/../Dropbox/SMART-R/Rda/sample_size.bootstrap.list.rda")
#save(sample_size.bootstrap.list.design.2, file = "~/../Dropbox/SMART-R/Rda/sample_size.bootstrap.list.design.2.rda")
