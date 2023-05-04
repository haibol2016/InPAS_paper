library("ggplot2")
library("patchwork")
library("reshape2")
library(ggridges)
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa)")

distance <- read.delim("001.InPAS.DaPars.pAPA.dist.diff.txt", header = TRUE)

distance_long <- melt(distance, id.vars = "Tx_ID",
                                variable.name = "Method", 
                                value.name = "Distance")

pdf(file = "Comparison of distance of predicted CPSs to nearest PAC polyA sites.pdf",
    width = 5, height = 3)
ggplot(distance_long, aes(x= Distance, fill = Method))+
  geom_histogram(bins =100) +
  scale_x_continuous(limits = c(-5000, 5000), 
                     n.breaks = 10,
                     expand = expansion(mult = 0, add = c(0, 100))) +
  scale_y_continuous(expand = expansion(mult = 0, add = c(0, 50))) +
  xlab("Distance to the nearest PAC-seq poly(A) sites (bp)") + 
  ylab("Count")+
  theme_classic() + theme(panel.border = 
                            element_rect(colour = "black", fill=NA, 
                                         linewidth = 0.1),
                          axis.title = element_text(size = 10),
                          axis.text = element_text(size = 8),
                          legend.key.size = unit(0.3, 'cm'))
dev.off()

