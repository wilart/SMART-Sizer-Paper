#Compute Power in the EXTEND Trial for IPW and AIPW
setwd("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/EXTEND")
V.IPW <- 250 * as.matrix(read.table("V.IPW.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = rep("numeric", 8))) #read in IPW matrix
V.DR <- 250 * as.matrix(read.table("V.DR.txt", sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = rep("numeric", 8))) #read in AIPW matrix

#write.table(round(V.IPW,2), file = "V.IPW.for.paper.txt", sep = "&")
#write.table(round(V.DR,2), file = "V.AIPW.for.paper.txt", sep = "&")

theta.IPW <- -c(7.56, 9.53, 8.05, 10.02, 7.71, 9.68, 8.19, 10.17)
theta.DR <- -c(7.65, 9.44, 7.83, 9.62, 8.06, 9.85, 8.24, 10.03)

Delta.IPW <- rep(max(theta.IPW), 8) - theta.IPW
Delta.DR <- rep(max(theta.DR), 8) - theta.DR

sapply(1:8, function(x) Delta.DR[x]/sqrt(V.DR[1,1]+V.DR[x,x]-2*V.DR[1,x]))


source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/simulate-normal.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-SVD.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-c-quantiles.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-power.R")
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/compute-sample-size.R")

set.seed(12478)
mc_list <- simulateNormal(n_EDTR = 8, n_sim = 1000)
mc_list_SVD_IPW <- ComputeSVD(V = V.IPW, mc_list)
mc_list_SVD_DR <- ComputeSVD(V = V.DR, mc_list)

c_alpha_IPW <- computeC(V = V.IPW, mc_list_SVD_IPW, alpha = 0.05)
c_alpha_DR <- computeC(V = V.DR, mc_list_SVD_DR, alpha = 0.05)


round(computePower(V.IPW, Delta = Delta.IPW, min_Delta = 2.15, sample_size = 250, c_alpha_IPW, Z = mc_list_SVD_IPW), 4) #compute IPW power
round(computePower(V.DR, Delta = Delta.DR, min_Delta = 2.15, sample_size = 250, c_alpha_DR, Z = mc_list_SVD_DR), 4) #compute AIPW power to eliminate all inferior DTRs

ceiling(computeSampleSize(V.IPW, Delta = Delta.IPW, min_Delta = 2.15, desired_power = 0.8, c_alpha_IPW, Z = mc_list_SVD_IPW))
ceiling(computeSampleSize(V.DR, Delta = Delta.DR, min_Delta = 2.15, desired_power = 0.8, c_alpha_DR, Z = mc_list_SVD_DR))

#############Plots

power.IPW.Delta.grid <- sapply(seq(0.14, 2.6, 0.01), function(x) computePower(V.IPW, Delta = Delta.IPW, min_Delta = x, sample_size = 250, c_alpha_IPW, Z = mc_list_SVD_IPW))
power.DR.Delta.grid <- sapply(seq(0.14, 2.6, 0.01), function(x) computePower(V.DR, Delta = Delta.DR, min_Delta = x, sample_size = 250, c_alpha_DR, Z = mc_list_SVD_DR))

library(ggplot2)
library(reshape2)

pdf(file = '~/../Dropbox/SMART-R/SMART-Power-Analysis/plots/EXTEND-power-by-delta-step3.pdf', width = 6, height = 5, onefile = FALSE)

power_melt_6 <- melt(data.frame(Delta = seq(0.14, 2.6, 0.01), IPW = power.IPW.Delta.grid, AIPW = power.DR.Delta.grid), id.vars = "Delta", variable.name = "Method")

g_power_6 <- ggplot(power_melt_6, aes(x = Delta, y = value, color =  Method, linetype = Method)) + geom_line() + theme_classic() +
  ggtitle(expression(paste("Power vs. ", Delta[min]))) +
  labs(x = expression(paste("Effect Size ", Delta[min])), y = "Power") + 
  theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size = 14))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2.6), breaks = seq(0,2.6,0.5)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.1))


g_power_6
dev.off()

power.IPW.equal.Delta.grid <- sapply(seq(0.14, 2.6, 0.01), function(x) computePower(V.IPW, Delta = c(0, rep(x,7)), min_Delta = x, sample_size = 250, c_alpha_IPW, Z = mc_list_SVD_IPW))
power.DR.equal.Delta.grid <- sapply(seq(0.14, 2.6, 0.01), function(x) computePower(V.DR, Delta = c(0, rep(x,7)), min_Delta = x, sample_size = 250, c_alpha_DR, Z = mc_list_SVD_DR))

power_melt_7 <- melt(data.frame(Delta = seq(0.14, 2.6, 0.01), IPW = power.IPW.equal.Delta.grid, AIPW = power.DR.equal.Delta.grid), id.vars = "Delta", variable.name = "Method")

pdf(file = '~/../Dropbox/SMART-R/SMART-Power-Analysis/plots/EXTEND-power-by-equal-delta5.pdf', width = 6, height = 5, onefile = FALSE)




g_power_equal_7 <- ggplot(power_melt_7, aes(x = Delta, y = value, color = Method, linetype = Method)) + geom_line() + theme_classic() +
  ggtitle("Power vs. Equal Effect Size") +
  labs(x = expression(paste("Effect Size ", Delta)), y = "Power") + 
  theme(axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 18), axis.text = element_text(size = 14))+
  scale_x_continuous(expand = c(0,0), limits = c(0, 2.6), breaks = seq(0,2.6,0.5)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.6), breaks = seq(0,0.6,0.1))



g_power_equal_7
dev.off()

#Compute actual sets of best
source("~/../Dropbox/SMART-R/SMART-Power-Analysis/R/Simulation/compute-best-set-empirical-power.R")
ComputeBestSet(theta = theta.IPW, c_alpha_IPW, V.IPW, sample_size = 250)
ComputeBestSet(theta = theta.DR, c_alpha_DR, V.DR, sample_size = 250)
