library(cowplot)
library(PASWR)
load("Figure5/Follicular_lymphoma/flfig.RData")
load("Figure5/Follicular_lymphoma/fl_bycell.RData")
load("Figure5/Normal_tonsil/tonsil.RData")
load("Figure5/GSE20300/GSE20300.RData")
load("Figure5/GSE20300/GSE20300_bycell.RData")


p1 = plot_grid(flfig, fl_bycell, ncol = 2, labels = c("A", "B"), label_size = 30)
p2 = plot_grid(tonsil, ncol = 1, labels = c("C"), label_size = 30)
p3 = plot_grid(GSE20300, gse20300.by.cell, ncol = 2, labels = c("D", "E"), label_size = 30)
plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1, 1, 1.2))