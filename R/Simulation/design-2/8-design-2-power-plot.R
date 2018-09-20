#generates the power plot for design 2, in figure 3
n_grid <- seq(100, 500, 50)

library(tidyr); library(dplyr); library(reshape2); library(RColorBrewer); library(ggplot2)

####
predicted.power <- melt(data.frame(n_grid, 'Predicted Power Identity' = design.2.predicted.power.identity,
                                   'Predicted Power True' = design.2.predicted.power.true, 
                                   'Predicted Power Conservative' = design.2.predicted.power.conservative,
                                   'Empirical Power Predicted' = design.2.empirical.power.by.n.delta.2,
                                   check.names = FALSE), id.vars = 'n_grid', value.name = 'Power')
levels(predicted.power$variable) <- c(expression(I[4]), expression(Sigma['True']), expression(Sigma['Conservative']), 'Empirical')



########

cl_predicted_power <- colorRampPalette(brewer.pal(n=6, name = 'Dark2'))(6)[c(1,2,4,3)]
g_power_by_n_predicted <- ggplot(predicted.power, aes(x = n_grid, y = Power, color = variable, linetype = variable)) + theme_classic() + 
  geom_point(show.legend = FALSE) + geom_line()  + xlab('Sample Size n') +ylab('Power') + 
  scale_color_manual(name = "", values = cl_predicted_power, guide = guide_legend(parse = TRUE), labels = c(expression(I[4]), expression(Sigma['True']), expression(Sigma['Conservative']), 'Empirical Power')) + 
  scale_linetype_manual(name = "", values = c(3,1,4,5), guide = guide_legend(parse = TRUE), labels = c(expression(I[4]), expression(Sigma['True']), expression(Sigma['Conservative']), 'Empirical Power')) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1.02), breaks = seq(0,1,0.1)) + 
  scale_x_continuous(expand = c(0,0), limits = c(75,510), breaks = seq(0,500,50)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.8) + theme(text = element_text(size = 14), 
                                       axis.title = element_text(size = 16), 
                                       axis.text = element_text(size = 14),
                                       legend.position = "bottom",
                                       legend.background = element_rect(fill = "white", size = 1, linetype = "solid", colour = "black")) + 
  ggtitle('Design 2: Power vs. Sample Size') +
  guides(colour = guide_legend(override.aes = list(shape = NA),keywidth = 1.4))

pdf(file = '~/../Dropbox/SMART-R/SMART-Power-Analysis/R/plots/SMART-design-2-predicted.pdf', width = 6, height = 5)
g_power_by_n_predicted
dev.off()



