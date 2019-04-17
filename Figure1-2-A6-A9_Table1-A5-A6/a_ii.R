library (nnls)
library (glmnet)
setwd("⁨Figure1-2_FigureA6-A9_Table1_TableA5-A6⁩")
source("sample_simulation.R")
source("../sourcecode/fardeep_function.R")
source("../sourcecode/Tuning_BIC.R")
source("../CIBERSORT/CIBERSORT.SVR.R")

compare = function(n = 500, p = 20, it = 50, a1 = 0.05, a2 = 0.2){
  iter = 1
  beta = matrix (0, it, p)
  coe.fardeep = matrix (0, it, p)
  coe.ciber = matrix (0, it, p)
  coe.nnls  = matrix (0, it, p)
  coe.dcq   = matrix (0, it, p)
  coe.ciber.remove = matrix (0, it, p)
  datax = NULL
  datay = matrix (0, n, it)
  tpr = NULL
  fpr = NULL
  para = NULL
  while(iter <= it) {
    sam = sample.sim(n = n, p = p, a1 = a1, a2 = a2, noise = "ii", outlier = "a", seed = iter + 50)
    x   = sam$x
    y   = sam$y
    datay[, iter] = y
    datax = cbind(datax, x)
    colnames(x) = c(1 : p)
    ### true beta
    beta[iter, ] = sam$beta
    ### cibersort
    cib.temp1 = CoreAlg(X = x, y = y)$w
    cib.id = which(cib.temp1 < 0)
    cib.temp1[cib.id] = 0
    coe.ciber[iter, ] = cib.temp1
    ### nnls
    coe.nnls [iter, ] = nnls(x, y)$x
    ### dcq
    elastic = glmnet(x, y, alpha = 0.05, lambda.min.ratio = 0.2, family = c('gaussian'), nlambda = 100,
                     intercept = FALSE)$beta
    dcq.temp1 = elastic[, ncol(elastic)]
    dcq.id = which(dcq.temp1 < 0)
    dcq.temp1[dcq.id] = 0
    coe.dcq[iter, ] = dcq.temp1
    ### fardeep
    kk = tuningBIC(x = x, y = y, n = 500, p = 20, intercept = TRUE, lognorm = FALSE)
    reg.fardeep = fardeep(x = x, y = y, k = kk, intercept = TRUE)
    coe.fardeep[iter, ] = reg.fardeep$beta[-1]
    loc  = sam$loc
    O    = a1 * n
    out  = sort (reg.fardeep$outlier_detect)
    para = c(para, kk)
    tpr  = c(tpr, 1 - sum (is.na(match (loc, out))) / O)
    fpr  = c(fpr, sum (is.na(match (out, loc))) / (n - O))
    ### ciber remove outlier
    newy = reg.fardeep$Y.new
    newx = reg.fardeep$X.new
    cibr = CoreAlg(X = newx, y = newy)$w
    cibr[cibr < 0] = 0
    coe.ciber.remove[iter, ] = cibr
    iter = iter + 1
  }
  tpr.fardeep = sum(tpr)/length(tpr)
  fpr.fardeep = sum(fpr)/length(fpr)
  sse.ciber = apply(coe.ciber - beta, 1, function(x) sum(x ^ 2))
  sse.fardeep = apply(coe.fardeep - beta, 1, function(x) sum(x ^ 2))
  sse.nnls = apply(coe.nnls - beta, 1, function(x) sum(x ^ 2))
  sse.dcq  = apply(coe.dcq  - beta, 1, function(x) sum(x ^ 2))
  sse.ciberr  = apply(coe.ciber.remove - beta, 1, function(x) sum(x ^ 2))
  write.table(datax, file= paste("sim_ii_a/datax", a1, ".txt", sep = ""), row.names=F, col.names=F, quote=F)
  write.table(datay, file= paste("sim_ii_a/datay", a1, ".txt", sep = ""), row.names=F, col.names=F, quote=F)
  write.table(beta, file= paste("sim_ii_a/beta", a1, ".txt", sep = ""), row.names=F, col.names=F, quote=F)
  result = list(tpr = tpr.fardeep, fpr = fpr.fardeep, k = para, sse.ciber = sse.ciber, sse.fardeep = sse.fardeep, sse.nnls = sse.nnls, sse.dcq = sse.dcq, sse.ciberr =sse.ciberr,
                mean.fardeep = mean(sse.fardeep), sd.fardeep = sd(sse.fardeep), mean.ciber = mean(sse.ciber),
                sd.ciber = sd(sse.ciber), mean.nnls = mean(sse.nnls), sd.nnls = sd(sse.nnls),
                mean.dcq = mean(sse.dcq), sd.dcq = sd(sse.dcq), coe.ciber = coe.ciber, coe.fardeep = coe.fardeep,
                coe.dcq = coe.dcq, coe.nnls = coe.nnls, coe.ciber.remove = coe.ciber.remove, coe = beta)
  return(result)
}

re_5  = compare(n = 500, p = 20, it = 50, a1 = 0.05, a2 = 0.2)
re_10 = compare(n = 500, p = 20, it = 50, a1 = 0.10, a2 = 0.2)
re_20 = compare(n = 500, p = 20, it = 50, a1 = 0.20, a2 = 0.2)
re_30 = compare(n = 500, p = 20, it = 50, a1 = 0.30, a2 = 0.2)

############ We got the reults of PERT by using the code from https://doi.org/10.1371/journal.pcbi.1002838.s003
sse.pert = read.table("../pert/sse_pert_ii_a.txt")
coe.pert = read.table("../pert/coe_pert_ii_a.txt")


## all methods
library(RColorBrewer)
library(wesanderson)
colo = c("grey56", "lightcoral", "seagreen4", "darkorchid4", "red")
require(ggplot2)
sse = matrix(c(re_5$sse.ciber, re_10$sse.ciber, re_20$sse.ciber, re_30$sse.ciber, re_5$sse.fardeep, 
               re_10$sse.fardeep, re_20$sse.fardeep, re_30$sse.fardeep, re_5$sse.nnls, re_10$sse.nnls, 
               re_20$sse.nnls, re_30$sse.nnls, re_5$sse.dcq, re_10$sse.dcq, re_20$sse.dcq,
               re_30$sse.dcq, sse.pert[, 1], sse.pert[, 2], sse.pert[, 3], sse.pert[, 4]), 1000, 1)
require(ggplot2)
dat = data.frame(Method = factor(rep(c("CIBERSORT","FARDEEP", "NNLS", "DCQ", "PERT"), each=200)),
                 Proportion.outlier = factor(rep(c("5% outlier", "10% outlier",
                                                   "20% outlier", "30% outlier", "5% outlier", "10% outlier",
                                                   "20% outlier", "30% outlier", "5% outlier", "10% outlier",
                                                   "20% outlier", "30% outlier", "5% outlier", "10% outlier",
                                                   "20% outlier", "30% outlier", "5% outlier", "10% outlier",
                                                   "20% outlier", "30% outlier"),each = 50)), SSE = log(sse, 2))
write.table (dat, file="sse.all.methods.log.txt", sep="\t", row.names=F, col.names=F, quote=F)
dat$Proportion.outlier = factor(dat$Proportion.outlier, levels=c("5% outlier", "10% outlier",
                                                                 "20% outlier", "30% outlier"))
dat$Method = factor(dat$Method, levels=c("NNLS", "DCQ", "PERT", "CIBERSORT", "FARDEEP"))
ii_a1 = ggplot(data = dat, aes(x = Proportion.outlier, y = SSE)) + geom_boxplot(aes(fill = Method), position=position_dodge(1)) + labs(x = "Proportion of outliers", y = "log2(SSE)") + scale_x_discrete(breaks = NULL) +
  facet_wrap( ~ Proportion.outlier, scales="free") + coord_cartesian(ylim=c(-7, 11)) + theme_gray() + scale_fill_manual(values=colo) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text=element_text(size=15), axis.title=element_text(size=15), strip.text = element_text(size=15))

save(ii_a1, file = "ii_a1.RData")

#######################

library(RColorBrewer)
library(gridExtra)
colo = c("grey56", "lightcoral", "seagreen4", "darkorchid4", "red")
dat = function(re = re, num = 1 : 50, per){
  b1 = as.vector(re$coe.fardeep)
  b2 = as.vector(re$coe.ciber)
  b3 = as.vector(re$coe.nnls)
  b4 = as.vector(re$coe.dcq)
  b5 = as.vector(as.matrix(coe.pert[num, ]))
  b6 = as.vector(re$coe)
  dd = cbind(b1, b2, b3, b4, b5)
  r2 = matrix(0, 1, 5)
  r = matrix(0, 1, 5)
  r2[1] = round(1 - sum((b1 - b6)^2)/sum((b6 - mean(b6))^2), 2)
  r2[2] = round(1 - sum((b2 - b6)^2)/sum((b6 - mean(b6))^2), 2)
  r2[3] = round(1 - sum((b3 - b6)^2)/sum((b6 - mean(b6))^2), 2)
  r2[4] = round(1 - sum((b4 - b6)^2)/sum((b6 - mean(b6))^2), 2)
  r2[5] = round(1 - sum((b5 - b6)^2)/sum((b6 - mean(b6))^2), 2)
  r[1] = round(cor(b1, b6), 2)
  r[2] = round(cor(b2, b6), 2)
  r[3] = round(cor(b3, b6), 2)
  r[4] = round(cor(b4, b6), 2)
  r[5] = round(cor(b5, b6), 2)
  meth = rbind(as.matrix(b1), as.matrix(b2), as.matrix(b3), as.matrix(b4), as.matrix(b5))
  true = as.matrix(rep(b6, 5))
  dd = data.frame(Method = factor(rep(c(paste("FARDEEP (R2 = ", r2[1],", R = ", r[1], ")"), 
                                        paste("CIBERSORT (R2 = ", r2[2],", R = ", r[2], ")"), 
                                        paste("NNLS (R2 = ", r2[3],", R = ", r[3], ")"),
                                        paste("DCQ (R2 = ", r2[4],", R = ", r[4], ")"), 
                                        paste("PERT (R2 = ", r2[5],", R = ", r[5], ")")), each=1000)), est = meth, tru = true)
  
  dd$Method = factor(dd$Method, levels= c(paste("NNLS (R2 = ", r2[3],", R = ", r[3], ")"), 
                                          paste("DCQ (R2 = ", r2[4],", R = ", r[4], ")"),
                                          paste("PERT (R2 = ", r2[5],", R = ", r[5], ")"),
                                          paste("CIBERSORT (R2 = ", r2[2],", R = ", r[2], ")"), 
                                          paste("FARDEEP (R2 = ", r2[1],", R = ", r[1], ")")))
  f = ggplot(dd, aes(est, tru, colour = Method, fill = Method)) + geom_point(cex = 0.5) +
    labs(x = "Estimated coefficients", y = "True coefficients") + ggtitle(per) +
    coord_cartesian(xlim=c(0, 1)) + guides(colour = guide_legend(override.aes = list(size=1))) + theme_gray() + scale_colour_manual(values = colo) +
    theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.position = c(0.75, 0.15), plot.title = element_text(hjust = 0.5), legend.background = element_rect(color = "black", fill = "grey90", size = 0.5, linetype = "solid"))
  return(list(dat = dd, fig = f))
}

res1 = dat(re = re_20, num = 101 : 150, per = "20% Outliers")
res2 = dat(re = re_30, num = 151 : 200, per = "30% Outliers")
res4 = dat(re = re_5, num = 1 : 50, per = "5% Outliers")
res5 = dat(re = re_10, num = 51 : 100, per = "10% Outliers")

library(cowplot)
library(PASWR)
ii_a2 = plot_grid(res4$fig, res5$fig, res1$fig, res2$fig, ncol = 2)

save(ii_a2, file = "ii_a2.RData")
