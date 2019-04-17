setwd("⁨Figure1-2_FigureA6-A9_Table1_TableA5-A6⁩")
library(cowplot)
library(PASWR)

load("i_a1.RData")
load("ii_a1.RData")

p1 = plot_grid(i_a1, ii_a1, ncol = 2, labels = c("A", "B"), label_size = 30)



load("i_a2.RData")
load("ii_a2.RData")


p2 = plot_grid(i_a2, ii_a2, ncol = 2, labels = c("A", "B"), label_size = 20)


load("i_b1.RData")
load("ii_b1.RData")

p3 = plot_grid(i_b1, ii_b1, ncol = 2, labels = c("A", "B"), label_size = 30)



load("i_b2.RData")
load("ii_b2.RData")


p4 = plot_grid(i_b2, ii_b2, ncol = 1, labels = c("A", "B"), label_size = 20)

