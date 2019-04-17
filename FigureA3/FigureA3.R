############# OV data
surv.ov = read.table("Data/OV_clin.txt", header = TRUE, row.names=1, sep = "\t", stringsAsFactors=F)
micro.ov = read.table("Data/OV.txt", sep = "\t", header = TRUE, row.names=1)
temp1 = colnames(micro.ov)
temp2 = colnames(surv.ov)
temp3 = strsplit(temp1, split = "[.]")
temp41 = lapply(temp3, function(x) x[1:3])
temp42 = as.vector(unlist(lapply(temp3, function(x) x[4])))
micro.tmp = micro.ov[, -c(grep("02", temp42), grep("11", temp42))]
temp4 = temp41[-c(grep("02", temp42), grep("11", temp42))]
temp5 = unlist(lapply(temp4, function(x) paste(x[1], x[2], x[3], sep = ".")))
temp6 = gsub("TCGA", "tcga", temp5)
sur.tmp.ov = surv.ov[, temp2 %in% temp6]
mic.tmp.ov = micro.tmp[, temp6 %in% temp2]
sur.ov = sur.tmp.ov[, order(colnames(sur.tmp.ov))]
mic.ov = mic.tmp.ov[, order(colnames(mic.tmp.ov))]

ge.ov = data.frame(gene = rownames(mic.ov), exp = 2 ^ mic.ov)
lm22 = read.table("Data/LM22.txt", header = TRUE, sep = "\t")
name.s = lm22$Gene.symbol
name.m = ge.ov$gene
gep.ov = ge.ov[name.m %in% name.s, ]
idm = order(gep.ov$gene)
gexp.ov = gep.ov[idm, -1]
row.names(gexp.ov) = gep.ov$gene[idm]
sg.ov = lm22[name.s %in% name.m, ]
ids = order(sg.ov$Gene.symbol)
sig.ov = sg.ov[ids, ]

################################## FARDEEP
source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
n     = nrow(gexp.ov)
p     = ncol(sig.ov) - 1
n.col = ncol(gexp.ov)
beta.fardeep.ov = matrix(0, n.col, p)
for (i in 1 : n.col){
  y = gexp.ov[, i]
  x = as.matrix(sig.ov[, -1])
  k  = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  reg   = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  beta.fardeep.ov[i, ] = reg$beta[-1]
  cat('i =', i, '\n')
}

purity = read.csv("Data/TCGA_purity.csv")
fardeep.immu.score.ov = as.data.frame(apply(beta.fardeep.ov, 1, sum))
type = sapply(strsplit(as.character(purity$Sample.ID), "-"), function(x) x[4])
purity = purity[type == "01A",]
samp1 = sapply(strsplit(as.character(purity$Sample.ID), "-"), function(x) paste0(c(x[1], x[2], x[3]), collapse = "."))
samp2 = sapply(strsplit(colnames(gexp.ov), "[.]"), function(x) paste0(c(x[2], x[3], x[4]), collapse = "."))
fardeep.immu.score.ov$Sample.ID = samp2
purity$Sample.ID = samp1
idx = samp1 %in% samp2
cpe.ov = purity[idx, c(1, 7)]
idx2 = samp2 %in% samp1
fardeep.score.ov = fardeep.immu.score.ov[idx2,]
na = which(is.na(cpe.ov$CPE))
r.ov = round(cor(fardeep.score.ov$`apply(beta.fardeep.ov, 1, sum)`[-na], cpe.ov$CPE[-na]), 2)
dat.ov = data.frame(FARDEEP.score = fardeep.score.ov$`apply(beta.fardeep.ov, 1, sum)`[-na], 
                    CPE = cpe.ov$CPE[-na], Source = rep(paste("OV (R =", r.ov, ")", sep = ""), length(cpe.ov$CPE[-na])))



############# LUSC data
surv.lusc = read.table("Data/LUSC_clin.txt", header = TRUE, row.names=1, sep = "\t", stringsAsFactors=F)
micro.lusc = read.table("Data/LUSC.txt", sep = "\t", header = TRUE, row.names=1)
temp1 = colnames(micro.lusc)
temp2 = colnames(surv.lusc)
temp3 = strsplit(temp1, split = "[.]")
temp4 = unlist(lapply(temp3, function(x) paste(x[1], x[2], x[3], sep = ".")))
temp5 = gsub("TCGA", "tcga", temp4)
sur.tmp.lusc = surv.lusc[, temp2 %in% temp5]
mic.tmp.lusc = micro.lusc[, temp5 %in% temp2]
sur.lusc = sur.tmp.lusc[, order(colnames(sur.tmp.lusc))]
mic.lusc = mic.tmp.lusc[, order(colnames(mic.tmp.lusc))]

ge.lusc = data.frame(gene = rownames(mic.lusc), exp = 2 ^ mic.lusc)
lm22 = read.table("Data/LM22.txt", header = TRUE, sep = "\t")
name.s = lm22$Gene.symbol
name.m = ge.lusc$gene
gep.lusc = ge.lusc[name.m %in% name.s, ]
idm = order(gep.lusc$gene)
gexp.lusc = gep.lusc[idm, -1]
row.names(gexp.lusc) = gep.lusc$gene[idm]
sg.lusc = lm22[name.s %in% name.m, ]
ids = order(sg.lusc$Gene.symbol)
sig.lusc = sg.lusc[ids, ]

################################## FARDEEP
source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
n     = nrow(gexp.lusc)
p     = ncol(sig.lusc) - 1
n.col = ncol(gexp.lusc)
beta.fardeep.lusc = matrix(0, n.col, p)
for (i in 1 : n.col){
  y = gexp.lusc[, i]
  x = as.matrix(sig.lusc[, -1])
  k  = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  reg   = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  beta.fardeep.lusc[i, ] = reg$beta[-1]
  cat('i =', i, '\n')
}

purity = read.csv("Data/TCGA_purity.csv")
fardeep.immu.score.lusc = as.data.frame(apply(beta.fardeep.lusc, 1, sum))
type = sapply(strsplit(as.character(purity$Sample.ID), "-"), function(x) x[4])
purity = purity[type == "01A",]
samp1 = sapply(strsplit(as.character(purity$Sample.ID), "-"), function(x) paste0(c(x[1], x[2], x[3]), collapse = "."))
samp2 = sapply(strsplit(colnames(gexp.lusc), "[.]"), function(x) paste0(c(x[2], x[3], x[4]), collapse = "."))
fardeep.immu.score.lusc$Sample.ID = samp2
purity$Sample.ID = samp1
idx = samp1 %in% samp2
cpe.lusc = purity[idx, c(1, 7)]
idx2 = samp2 %in% samp1
fardeep.score.lusc = fardeep.immu.score.lusc[idx2,]
r.lusc = round(cor(fardeep.score.lusc$`apply(beta.fardeep.lusc, 1, sum)`, cpe.lusc$CPE), 2)
dat.lusc = data.frame(FARDEEP.score = fardeep.score.lusc$`apply(beta.fardeep.lusc, 1, sum)`, 
                      CPE = cpe.lusc$CPE, Source = rep(paste("LUSC (R =", r.lusc, ")", sep = ""), length(cpe.lusc$CPE)))
library(ggplot2)
library(grid)
dat = rbind(dat.ov, dat.lusc)
p = ggplot(dat, aes(CPE, FARDEEP.score)) +
  geom_point(cex = 1, color = "blue") + labs(x = "CPE", y = "FARDEEP SCORE") + geom_smooth(method='lm', color = "black", se=FALSE) + 
  facet_wrap( ~ Source, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))
p + theme(panel.spacing = unit(2, "lines"))

