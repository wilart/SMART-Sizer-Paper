#Simulation Design 1: generates the power plot

n_grid <- seq(100, 500, 50)

library(reshape2); library(ggplot2); library(RColorBrewer); library(dplyr)

design.1.predicted.power <- melt(data.frame(n_grid,'Predicted Power Identity' = design.1.predicted.power.identity, 'Predicted Power True' = design.1.predicted.power.true, 'Predicted Power No Correlation' = design.1.predicted.power.nocor, 'True Power' = design.1.empirical.power.by.n.delta.0.25, check.names = FALSE), id.vars = 'n_grid', value.name = 'Power')

levels(design.1.predicted.power$variable) <- c(expression(I[4]), expression(Sigma['True']), expression(paste(sigma['max']^2,I['4'])), 'True/Power')

cl_predicted_power <- colorRampPalette(brewer.pal(n=6, name = 'Dark2'))(6)[c(1,2,4,3)]
g_power_by_n_predicted <- ggplot(design.1.predicted.power, aes(x = n_grid, y = Power, color = variable, linetype = variable)) + 
  theme_classic() + 
  geom_point() + 
  geom_line()  + 
  xlab('Sample Size n') +
  ylab('Power') + 
  scale_color_manual(name = "", values = cl_predicted_power, guide = guide_legend(parse = TRUE), labels = c(expression(I[4]), expression(Sigma['True']),expression(Sigma['Conservative']),'Empirical Power')) + 
  scale_linetype_manual(name = "", values = c(3,1,4,5), guide = guide_legend(parse = TRUE), labels = c(expression(I[4]), expression(Sigma['True']),expression(Sigma['Conservative']),'Empirical Power')) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.02), breaks = seq(0,1,0.1)) + 
  scale_x_continuous(expand = c(0,0), limits = c(75,510), breaks = seq(0,500,50)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.8) + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", colour = "black"),
        legend.title =element_text(size = 14),
        legend.position = "bottom") + 
  ggtitle("Design 1: Power vs. Sample Size") + guides(colour = guide_legend(override.aes = list(shape = NA)))+
  guides(colour = guide_legend(override.aes = list(shape = NA),keywidth = 1.4))



pdf(file = '~/../Dropbox/SMART-R/SMART-Power-Analysis/R/Plots/SMART-design-1-predicted.pdf', width = 6, height = 5)

g_power_by_n_predicted

dev.off()
