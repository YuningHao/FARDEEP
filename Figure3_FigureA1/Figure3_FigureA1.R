# library(annotate)
# library(preprocessCore)
# library(AnnotationDbi)
# library(affy)
# source("sourcecode/fardeep_function.R")
# source("sourcecode/Tuning_BIC.R")
# source("CIBERSORT/CIBERSORT_no_last_normalization.R")
# ee = NULL
# for (i in 279601 : 279612){
#   cel = ReadAffy(filenames = paste("Data/GSE11058_RAW/GSM", i, ".CEL.gz", sep = ""))
#   e   = exprs(mas5(cel))
#   ee = cbind(ee, e)
# }
# exp = normalize.quantiles(ee)
# ex = data.frame(probe = rownames(ee), ex = exp)
# write.csv(ex, "Figure3_FigureA1/GSE11103_mix.csv", row.names = FALSE)


source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
source("CIBERSORT/CIBERSORT_no_last_normalization.R")
ex = read.csv("Figure3_FigureA1/GSE11103_mix.csv")
sig = read.table("Data/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt",
                 header = TRUE) 
name = data.frame(probe = sig$NAME)
gep = merge(name, ex, by = "probe")[ , -1]
row.names(gep) = name$probe

mixA   = c(0.25, 0.125, 0.25, 0.375)
mixB   = c(0.05, 0.317, 0.475, 0.158)
mixC   = c(0.01, 0.495, 0.165, 0.33)
mixD   = c(0.002, 0.333, 0.333, 0.333)
true.beta = t (matrix(c(rep(mixA, 3), rep(mixB, 3), rep(mixC, 3), rep(mixD, 3)), 4, 12))


m = gep
s = sig[, -1]
rownames(s) = sig[, 1]
######### FARDEEP
n     = nrow(m)
p     = ncol(s)
n.col = ncol(m)
beta.fardeep = matrix(0, n.col, p)
para  = NULL
nout = NULL
outlier = matrix(0, n, n.col)
for (i in 1 : 12){
  y = m [, i]
  x = as.matrix(s)
  k = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  para    = rbind (para, k)
  reg     = fardeep(x = x, y = y, k = k, intercept = TRUE)
  nout    = c(nout, reg$number_outlier)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  coe     = reg$beta[-1]
  beta.fardeep[i, ] = coe
}
diff.fardeep = abs(beta.fardeep - true.beta)
sse.fardeep  = apply(diff.fardeep, 1, function(x) sum (x^2))

#out_gene = NULL
#for (i in 1:ncol(outlier)){
#  out_tmp = sig$NAME[outlier[, i]]
#  out_gene = c(out_gene, paste0(out_tmp, collapse = ", "))
#}
#out_dat = data.frame(Samples = paste("GSM", 279601 : 279612, sep = ""), Number.of.outlier = nout, Outlier.probes = out_gene)
#write.csv(out_dat, "number_outlier_abbas.csv", row.names = FALSE)

######### CIBERSORT
write.table(ex, "mixture_file.txt", sep="\t", row.names = FALSE)
write.table(sig, "sig_matrix.txt", sep="\t", row.names = FALSE)

result = CIBERSORT.n(sig_matrix = "sig_matrix.txt" , 
                   mixture_file = "mixture_file.txt", perm=0, QN=F)
beta.ciber = result[ , 1:4]
diff.ciber = abs(beta.ciber - true.beta)
sse.ciber  = apply (diff.ciber, 1, function(x) sum (x^2))


########### CIBERSORT absolute
beta.ciberabs = read.table ("Data/CIBERSORT.abs.txt", sep = "\t", header = TRUE)[, 1:4]
diff.ciberabs = abs(beta.ciberabs - true.beta)
sse.ciberabs  = apply(diff.ciberabs, 1, function(x) sum(x^2))


########## NNLS
n     = nrow(m)
p     = ncol(s)
n.col = ncol(m)
beta.nnls = matrix(0, n.col, p)
para  = NULL
for (i in 1 : n.col){
  y = m [, i]
  x = as.matrix(s)
  model = nnls (x, y)
  coe   = model$x
  beta.nnls[i, ] = coe
}
diff.nnls = abs(beta.nnls - true.beta)
sse.nnls  = apply (diff.nnls, 1, function(x) sum (x^2))


########## DCQ
library(glmnet)
n     = nrow(m)
p     = ncol(s)
n.col = ncol(m)
beta.dcq = matrix(0, n.col, p)
para  = NULL
for (i in 1 : n.col){
  y = m[, i]
  x = as.matrix(s)
  model = glmnet (x, y, alpha = 0.05, lambda.min.ratio = 0.2, family = c('gaussian'), nlambda = 100,
                  intercept = FALSE)$beta
  coe   = model[, ncol(model)]
  beta.dcq[i, ] = coe
}
diff.dcq = abs(beta.dcq - true.beta)
sse.dcq  = apply (diff.dcq, 1, function(x) sum (x^2))



########## PERT
#write.table(gep, "./perty.txt", sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(s, "./pertx.txt", sep="\t", row.names = FALSE, col.names = FALSE)
beta.pert = read.table ("pert/pert_figure3_figureA1.txt")
diff.pert = abs(beta.pert - true.beta)
sse.pert  = apply (diff.pert, 1, function(x) sum (x^2))

library(ggplot2)
library(plotly)
library(wesanderson)
library(RColorBrewer)
colo = c("grey56", "lightcoral", "seagreen4", "darkorchid4", "red", "blue")
sse = t(cbind(t(sse.fardeep), t(sse.ciber), t(sse.pert), t(sse.nnls), t(sse.dcq)))

dat = data.frame(Method = factor(rep(c("FARDEEP", "CIBERSORT", "PERT", "NNLS", "DCQ"), each=12)), 
                 MIX = factor(rep(c("mixA1", "mixA2","mixA3","mixB1","mixB2","mixB3","mixC1","mixC2","mixC3",
                                    "mixD1","mixD2","mixD3"), 5)), cat = factor(rep(c("mixA", "mixA","mixA","mixB",
                                                                                      "mixB","mixB","mixC","mixC","mixC","mixD","mixD","mixD"), 5)), SSE = sse)
dat$Method = factor(dat$Method, levels=c("NNLS", "DCQ", "PERT", "CIBERSORT", "FARDEEP"))

p1 = ggplot(data=dat, aes(x = MIX, y = SSE, fill = Method)) + labs(x = "") + 
  geom_bar(colour = "black", size = 0.2, stat = "identity", position=position_dodge()) + ylim(low = 0, high = 0.22) + facet_wrap( ~ cat, scales="free") + scale_fill_manual(values=colo) +
  theme_gray() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 15), axis.title=element_text(size=15), axis.text=element_text(size=14), strip.text = element_text(size=15))


############## r2
af = as.vector(beta.fardeep)
ac = as.vector(beta.ciber)
an = as.vector(beta.nnls)
ad = as.vector(beta.dcq)
ap = as.vector(as.matrix(beta.pert))
at = as.vector(true.beta)

r2.fardeep = round(1 - sum((af-at)^2)/sum((at - mean(at))^2), 2)
r2.ciber = round(1 - sum((ac-at)^2)/sum((at - mean(at))^2), 2)
r2.nnls = round(1 - sum((an-at)^2)/sum((at - mean(at))^2), 2)
r2.pert = round(1 - sum((ap-at)^2)/sum((at - mean(at))^2), 2)
r2.dcq = round(1 - sum((ad-at)^2)/sum((at - mean(at))^2), 2)

r.fardeep = round(cor(af, at), 2)
r.ciber = round(cor(ac, at), 2)
r.nnls = round(cor(an, at), 2)
r.pert = round(cor(ap, at), 2)
r.dcq = round(cor(ad, at), 2)

coe.tru = rep(at, 5)
coe.est = c(af, ac, an, ap, ad)
method  = rep(c(paste("FARDEEP (R2 = ", r2.fardeep, ", R = ", r.fardeep, ")", sep = ""), 
                paste("CIBERSORT (R2 = ", r2.ciber, ", R = ", r.ciber, ")", sep = ""), 
                paste("NNLS (R2 = ", r2.nnls, ", R = ", r.nnls, ")", sep = ""),
                paste("PERT (R2 = ", r2.pert, ", R = ", r.pert, ")", sep = ""),
                paste("DCQ (R2 = ", r2.dcq, ", R = ", r.dcq, ")", sep = "")), each = 48)
cell.type = rep(rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 12), 5)
dat = data.frame(Method = factor(method), Cell.type = cell.type, coe_est = coe.est, coe_tru = coe.tru)
dat$Method = factor(dat$Method, levels = c(paste("NNLS (R2 = ", r2.nnls, ", R = ", r.nnls, ")", sep = ""), 
                                           paste("DCQ (R2 = ", r2.dcq, ", R = ", r.dcq, ")", sep = ""),
                                           paste("PERT (R2 = ", r2.pert, ", R = ", r.pert, ")", sep = ""), 
                                           paste("CIBERSORT (R2 = ", r2.ciber, ", R = ", r.ciber, ")", sep = ""), 
                                           paste("FARDEEP (R2 = ", r2.fardeep, ", R = ", r.fardeep, ")", sep = "")))


library(ggplot2)
p2  = ggplot(dat, aes(coe_est, coe_tru, shape = Cell.type, colour = Cell.type)) + xlim(0, 0.7) + ylim(0, 0.7) +
  geom_point(cex = 1.5) + labs(y = "True cell amount", x = "Estimated cell amount") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Method, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "right", legend.text = element_text(size = 15), axis.title=element_text(size=15), strip.text = element_text(size=12))


######### by cell
cf = round(c(cor(beta.fardeep[, 1], true.beta[, 1]), cor(beta.fardeep[, 2], true.beta[, 2]), 
             cor(beta.fardeep[, 3], true.beta[, 3]), cor(beta.fardeep[, 4], true.beta[, 4])), 2)
cc = round(c(cor(beta.ciber[, 1], true.beta[, 1]), cor(beta.ciber[, 2], true.beta[, 2]), 
             cor(beta.ciber[, 3], true.beta[, 3]), cor(beta.ciber[, 4], true.beta[, 4])), 2)
sf = round(c(sum((beta.fardeep[, 1] - true.beta[, 1])^2), sum((beta.fardeep[, 2] - true.beta[, 2])^2), 
             sum((beta.fardeep[, 3] - true.beta[, 3])^2), sum((beta.fardeep[, 4] - true.beta[, 4])^2)), 2)
sc = round(c(sum((beta.ciber[, 1] - true.beta[, 1])^2), sum((beta.ciber[, 2] - true.beta[, 2])^2), 
             sum((beta.ciber[, 3] - true.beta[, 3])^2), sum((beta.ciber[, 4] - true.beta[, 4])^2)), 2)
r2f = round(c(1 - sum((beta.fardeep[, 1] - true.beta[, 1])^2)/sum((true.beta[, 1] - mean(true.beta[, 1]))^2), 1 - sum((beta.fardeep[, 2] - true.beta[, 2])^2)/sum((true.beta[, 2] - mean(true.beta[, 2]))^2), 
              1 - sum((beta.fardeep[, 3] - true.beta[, 3])^2)/sum((true.beta[, 3] - mean(true.beta[, 3]))^2), 1 - sum((beta.fardeep[, 4] - true.beta[, 4])^2)/sum((true.beta[, 4] - mean(true.beta[, 4]))^2)), 2)
r2c = round(c(1 - sum((beta.ciber[, 1] - true.beta[, 1])^2)/sum((true.beta[, 1] - mean(true.beta[, 1]))^2), 1 - sum((beta.ciber[, 2] - true.beta[, 2])^2)/sum((true.beta[, 2] - mean(true.beta[, 2]))^2), 
              1 - sum((beta.ciber[, 3] - true.beta[, 3])^2)/sum((true.beta[, 3] - mean(true.beta[, 3]))^2), 1 - sum((beta.ciber[, 4] - true.beta[, 4])^2)/sum((true.beta[, 4] - mean(true.beta[, 4]))^2)), 2)

val = c(af, ac)
tru = rep(at, 2)
meth = rep(c("FARDEEP", "CIBERSORT"), each = 48)
cell = rep(rep(c("Jurkat", "IM.9", "Raji", "THP.1"), each = 12), 2)
dat.r = data.frame(Method = meth, Cell.type = cell, coe_est = val, coe_tru = tru)
library(ggplot2)
library(grid)
colo2 = c("darkorchid4", "red")
p3 = ggplot(dat.r, aes(coe_est, coe_tru, colour = Method)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(y = "True cell amount", x = "Estimated cell amount") + geom_abline(intercept = 0, slope = 1) + scale_colour_manual(values=colo2) + 
  facet_wrap( ~ Cell.type, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))



text1 = paste("FARDEEP \n (R = ", cf[1], ", SSE = ", sf[1], ", R2 = ", r2f[1], ") \n  CIBERSORT \n (R = ", cc[1], ", SSE = ", sc[1], ", R2 = ", r2c[1], ")", sep = "")
text2 = paste("FARDEEP \n (R = ", cf[2], ", SSE = ", sf[2], ", R2 = ", r2f[2], ") \n  CIBERSORT \n (R = ", cc[2], ", SSE = ", sc[2], ", R2 = ", r2c[2], ")", sep = "")
text3 = paste("FARDEEP \n (R = ", cf[3], ", SSE = ", sf[3], ", R2 = ", r2f[3], ") \n  CIBERSORT \n (R = ", cc[3], ", SSE = ", sc[3], ", R2 = ", r2c[3], ")", sep = "")
text4 = paste("FARDEEP \n (R = ", cf[4], ", SSE = ", sf[4], ", R2 = ", r2f[4], ") \n  CIBERSORT \n (R = ", cc[4], ", SSE = ", sc[4], ", R2 = ", r2c[4], ")", sep = "")

text = c(text1, text2, text3, text4)
cell = c("Jurkat", "IM.9", "Raji", "THP.1")
graphLabels = data.frame(Cell.type = cell, text = text)

by.cell = p3 + geom_text(data = graphLabels, aes(x = 0.25, y = 0.8, label = text, color = NULL), size = 3)



library(cowplot)
library(PASWR)
p4 = plot_grid(p1, p2, by.cell, ncol = 1, labels = c("A", "B", "C"), label_size = 20)




########## bar plot with abs
library(ggplot2)
library(plotly)
library(wesanderson)
library(RColorBrewer)
colo = c("grey56", "lightcoral", "seagreen4", "darkorchid4", "red", "blue")
yticks = c(0, 0.05, 0.1, 0.15, 0.2, 1, 2, 3, 4)
trans = function(x){pmin(x,0.22) + 0.05*pmax(x-0.22,0)}

sse = trans(t(cbind (t(sse.fardeep), t(sse.ciber), t(sse.ciberabs), t(sse.pert), t(sse.nnls), t(sse.dcq))))

dat = data.frame(Method = factor(rep(c("FARDEEP", "CIBERSORT", "CIBERSORT.abs", "PERT", "NNLS", "DCQ"), each=12)), 
                 MIX = factor(rep(c("mixA1", "mixA2","mixA3","mixB1","mixB2","mixB3","mixC1","mixC2","mixC3",
                                    "mixD1","mixD2","mixD3"), 6)), cat = factor(rep(c("mixA", "mixA","mixA","mixB",
                                                                                      "mixB","mixB","mixC","mixC","mixC","mixD","mixD","mixD"), 6)), SSE = sse)
dat$Method = factor(dat$Method, levels=c("NNLS", "DCQ", "PERT", "CIBERSORT", "CIBERSORT.abs",  "FARDEEP"))

pp1 = ggplot(data=dat, aes(x = MIX, y = SSE, fill = Method)) + labs(x = "") + 
  geom_bar(colour = "black", size = 0.2, stat = "identity", position=position_dodge()) + scale_colour_manual(values=colo) + 
  facet_wrap( ~ cat, scales="free") + scale_fill_manual(values=colo) + 
  theme_gray() + geom_rect(aes(xmin=0, xmax=4, ymin=0.22, ymax = 0.25), fill="white") + scale_y_continuous(limits=c(0,NA), breaks=trans(yticks), labels=yticks) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 15), axis.title=element_text(size=15), axis.text=element_text(size=14), strip.text = element_text(size=15))

write.csv(dat, "SSE_cell_line.csv")

############## r2
af = as.vector(beta.fardeep)
ac = as.vector(beta.ciber)
aca = as.vector(as.matrix(beta.ciberabs))
an = as.vector(beta.nnls)
ad = as.vector(beta.dcq)
ap = as.vector(as.matrix(beta.pert))
at = as.vector(true.beta)
r2.fardeep = round(1 - sum((af-at)^2)/sum((at - mean(at))^2), 2)
r2.ciber = round(1 - sum((ac-at)^2)/sum((at - mean(at))^2), 2)
r2.ciberabs = round(1 - sum((aca-at)^2)/sum((at - mean(at))^2), 2)
r2.nnls = round(1 - sum((an-at)^2)/sum((at - mean(at))^2), 2)
r2.pert = round(1 - sum((ap-at)^2)/sum((at - mean(at))^2), 2)
r2.dcq = round(1 - sum((ad-at)^2)/sum((at - mean(at))^2), 2)

r.fardeep = round(cor(af, at), 2)
r.ciber = round(cor(ac, at), 2)
r.ciberabs = round(cor(aca, at), 2)
r.nnls = round(cor(an, at), 2)
r.pert = round(cor(ap, at), 2)
r.dcq = round(cor(ad, at), 2)

coe.tru = rep(at, 6)
coe.est = c(af, ac, aca, an, ap, ad)
method  = rep(c(paste("FARDEEP (R2 = ", r2.fardeep, ", R = ", r.fardeep, ")", sep = ""), 
                paste("CIBERSORT (R2 = ", r2.ciber, ", R = ", r.ciber, ")", sep = ""), 
                paste("CIBERSORT.abs (R2 = ", r2.ciberabs, ", R = ", r.ciberabs, ")", sep = ""),
                paste("NNLS (R2 = ", r2.nnls, ", R = ", r.nnls, ")", sep = ""),
                paste("PERT (R2 = ", r2.pert, ", R = ", r.pert, ")", sep = ""),
                paste("DCQ (R2 = ", r2.dcq, ", R = ", r.dcq, ")", sep = "")), each = 48)
cell.type = rep(rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 12), 6)
dat = data.frame(Method = factor(method), Cell.type = cell.type, coe_est = coe.est, coe_tru = coe.tru)
dat$Method = factor(dat$Method, levels = c(paste("NNLS (R2 = ", r2.nnls, ", R = ", r.nnls, ")", sep = ""), 
                                           paste("DCQ (R2 = ", r2.dcq, ", R = ", r.dcq, ")", sep = ""),
                                           paste("PERT (R2 = ", r2.pert, ", R = ", r.pert, ")", sep = ""), 
                                           paste("CIBERSORT (R2 = ", r2.ciber, ", R = ", r.ciber, ")", sep = ""),
                                           paste("CIBERSORT.abs (R2 = ", r2.ciberabs, ", R = ", r.ciberabs, ")", sep = ""),
                                           paste("FARDEEP (R2 = ", r2.fardeep, ", R = ", r.fardeep, ")", sep = "")))
library(ggplot2)
pp2  = ggplot(dat, aes(coe_est, coe_tru, shape = Cell.type, colour = Cell.type)) + ylim(0, 0.7) + xlim(0, 1.95) +
  geom_point(cex = 1.5) + labs(y = "True cell amount", x = "Estimated cell amount") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Method, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.text = element_text(size = 15), axis.title=element_text(size=15), strip.text = element_text(size=12))

library(cowplot)
library(PASWR)

pp3 = plot_grid(pp1, pp2, ncol = 1, labels = c("A", "B"), label_size = 20)
