#Computes the power for all pairwise comparisons (Ogbagaber, 2016) and for MCB.
#To produce Table 1 in the supplementary materials

source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")

theta.DR <- c(1.802, 1.3, 1.699, 1.197)
Delta.DR <- max(theta.DR) - theta.DR
library(Matrix)
load("~/../Dropbox/SMART-R/Rda/design.1.cov.true.list.rda")
design.1.cov.true <- as.matrix(nearPD(apply(simplify2array(design.1.cov.true.list),1:2, mean))$mat)

V2 <- design.1.cov.true
V4 <- diag(diag(V2), 4)

sample_size_4 <- c(max(apply(combn(c(1,2,3,4), m = 2),2,function(x) (V2[x[1],x[1]]+V2[x[2],x[2]]-2*V2[x[1],x[2]])*(qnorm(1-0.05/(4*3))+qnorm(1-0.2))^2/(theta.DR[x[1]]-theta.DR[x[2]])^2)),
max(apply(combn(c(1,2,3,4), m = 2),2,function(x) (V4[x[1],x[1]]+V4[x[2],x[2]]-2*V4[x[1],x[2]])*(qnorm(1-0.05/(4*3))+qnorm(1-0.2))^2/(theta.DR[x[1]]-theta.DR[x[2]])^2)))

set.seed(17263)
mc_list <- simulateNormal(n_EDTR = 4, n_sim = 1000)
MCB_sample_size_4 <- 
  c(computeSampleSize(V = V2, Delta = Delta.DR, min_Delta = 0.1, c_alpha = computeC(V = V2, mc_list_SVD = ComputeSVD(V = V2, mc_list)), desired_power = 0.8, Z = ComputeSVD(V = V2, mc_list)),
    computeSampleSize(V = V4, Delta = Delta.DR, min_Delta = 0.1, c_alpha = computeC(V = V4, mc_list_SVD = ComputeSVD(V = V4, mc_list)), desired_power = 0.8, Z = ComputeSVD(V = V4, mc_list)))


###############Design 2
theta.DR <- c(1.499627, 3.500691, 3.251208, 4.251172, 3.50103)
Delta.DR <- max(theta.DR) - theta.DR

load('~/../Dropbox/SMART-R/Rda/design.2.cov.true.list.rda')
design.2.cov.true <- apply(simplify2array(design.2.cov.true.list),1:2, mean)
V2 <- design.2.cov.true
V4 <- diag(diag(V2), 5)

comparisons.all <- combn(c(1,2,3,4,5), m = 2)

rbind(1:ncol(comparisons.all),apply(comparisons.all,2,function(x) theta.DR[x[1]]-theta.DR[x[2]])) #consider pairwise comparisons which have true difference 0.7 or more
comparisons <- comparisons.all[,c(1,2,3,4,6,8,10)]


sample_size_5 <- c(max(apply(comparisons,2,function(x) (V2[x[1],x[1]]+V2[x[2],x[2]]-2*V2[x[1],x[2]])*(qnorm(1-0.05/(ncol(comparisons)*2))+qnorm(1-0.2))^2/(theta.DR[x[1]]-theta.DR[x[2]])^2)),
max(apply(comparisons,2,function(x) (V4[x[1],x[1]]+V4[x[2],x[2]]-2*V4[x[1],x[2]])*(qnorm(1-0.05/(ncol(comparisons)*2))+qnorm(1-0.2))^2/(theta.DR[x[1]]-theta.DR[x[2]])^2)))

set.seed(17263)
mc_list <- simulateNormal(n_EDTR = 5, n_sim = 1000)
MCB_sample_size_5 <- c(computeSampleSize(V = V2, Delta = max(theta.DR)-theta.DR, min_Delta = 0.7, c_alpha = computeC(V = V2, mc_list_SVD = ComputeSVD(V = V2, mc_list)), desired_power = 0.8, Z = ComputeSVD(V = V2, mc_list)),
computeSampleSize(V = V4, Delta = max(theta.DR)-theta.DR, min_Delta = 0.7, c_alpha = computeC(V = V4, mc_list_SVD = ComputeSVD(V = V4, mc_list)), desired_power = 0.8, Z = ComputeSVD(V = V4, mc_list)))


summary_table <- ceiling(cbind(sample_size_4, sample_size_5, MCB_sample_size_4, MCB_sample_size_5))
colnames(summary_table) <- c("4 Pairwise", "5 Pairwise", "4 MCB", "5 MCB")
rownames(summary_table) <- c("True Matrix", "Conservative")
library(knitr)
kable(summary_table,format = 'latex', caption = "Sample size to achieve power for pairwise comparisons and for MCB")

