library(cowplot)
library(PASWR)
load("Figure6_FigureS2/LUSC/lusc.RData")
load("Figure6_FigureS2/OV/OV.RData")

pcomb = plot_grid(ov, lusc, ncol = 1, labels = c("A", "B"), label_size = 20, rel_heights = c(1, 1))


load("Figure6_FigureS2/LUSC/LUSC_abs.RData")
load("Figure6_FigureS2/OV/OV_abs.RData")

pcomb2 = plot_grid(ov_abs, lusc_abs, ncol = 2, labels = c("A", "B"), label_size = 20, rel_heights = c(1, 1))
