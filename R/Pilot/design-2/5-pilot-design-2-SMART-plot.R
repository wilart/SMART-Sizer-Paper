#Simulation Design 2: generates the power plot for design 2

n_grid <- seq(100, 600, 50)

library(reshape2); library(ggplot2); library(RColorBrewer); library(dplyr)

pilot.power.melt <- data.frame(cbind(melt(data.frame(pilot.power.summary), id.vars = "n_grid", value.name = "Power"), rbind(cbind(rep(NA,11),rep(NA,11)),pilot.unknown.var.CI[,2:3],pilot.known.var.CI[,2:3],cbind(rep(NA,11),rep(NA,11)))))



cl_predicted_power <- colorRampPalette(brewer.pal(n=4, name = 'Set1'))(4)[c(4,1,2,3)]

g_power_by_n_predicted <- ggplot(pilot.power.melt, aes(x = n_grid, y = Power, color = variable, linetype = variable)) + 
  theme_classic() + 
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymin = X2.5.,ymax=X97.5.), width = 14)+ 
  xlab('Sample Size n') +
  ylab('Power') + 
  scale_color_manual(name = "", values = cl_predicted_power, guide = guide_legend(parse = TRUE), labels = c(expression(Sigma["True"]),expression(paste(Sigma["Pilot, Unknown Var."])),expression(Sigma["Pilot, Known Var."]),expression(Sigma["Conservative"]))) + 
  scale_linetype_manual(name = "", values = c(2,1,1,4), guide = guide_legend(parse = TRUE), labels = c(expression(Sigma["True"]),expression(paste(Sigma["Pilot, Unknown Var."])),expression(Sigma["Pilot, Known Var."]),expression(Sigma["Conservative"]))) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.02), breaks = seq(0,1,0.1)) + 
  scale_x_continuous(expand = c(0,0), limits = c(75,610), breaks = seq(0,600,50)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.8) + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", colour = "black"),
        legend.title =element_text(size = 14),
        legend.position = "bottom") + 
  ggtitle("Pilot Design 2: Power vs. Sample Size") + guides(colour = guide_legend(override.aes = list(shape = NA)))+
  guides(colour = guide_legend(override.aes = list(shape = NA),keywidth = 2))



#pdf(file = '~/../Dropbox/SMART-R/SMART-Power-Analysis/R/Plots/SMART-Pilot-Power2.pdf', width = 6, height = 5)

g_power_by_n_predicted

dev.off()
