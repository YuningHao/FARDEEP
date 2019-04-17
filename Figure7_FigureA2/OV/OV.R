surv = read.table("Data/OV_clin.txt", header = TRUE, row.names=1, sep = "\t", stringsAsFactors=F)
micro = read.table("Data/OV.txt", sep = "\t", header = TRUE, row.names=1)

##################################################
##### combine sample clinical data with microarray data  #######
##################################################
temp1 = colnames(micro)
temp2 = colnames(surv)
temp3 = strsplit(temp1, split = "[.]")
temp41 = lapply(temp3, function(x) x[1:3])
temp42 = as.vector(unlist(lapply(temp3, function(x) x[4])))
micro.tmp = micro[, -c(grep("02", temp42), grep("11", temp42))]
temp4 = temp41[-c(grep("02", temp42), grep("11", temp42))]
temp5 = unlist(lapply(temp4, function(x) paste(x[1], x[2], x[3], sep = ".")))
temp6 = gsub("TCGA", "tcga", temp5)
surr  = surv[, temp2 %in% temp6]
micc = micro.tmp[, temp6 %in% temp2]
sur = surr[, order(colnames(surr))]
mic = micc[, order(colnames(micc))]
write.table(mic, "./ovmerge.txt", sep="\t", row.names=TRUE, col.names = TRUE)

ge = data.frame(gene = rownames(mic), exp = 2 ^ mic)  # the first colum of ge is gene name
lm22 = read.table("Data/LM22.txt", header = TRUE, sep = "\t")
name.s = lm22$Gene.symbol
name.m = ge$gene
gep = ge[name.m %in% name.s, ]
idm = order(gep$gene)
gexp = gep[idm, -1]
row.names(gexp) = gep$gene[idm]
sg = lm22[name.s %in% name.m, ]
ids = order(sg$Gene.symbol)
sig = sg[ids, ]

################################## FARDEEP
source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
n     = nrow(gexp)
p     = ncol(sig) - 1
n.col = ncol(gexp)
beta.fardeep = matrix(0, n.col, p)
outlier = matrix(0, n, n.col)
nout  = NULL
para  = NULL
inter = 1:dim(gexp)[2]
for (i in 1 : n.col){
  y = gexp[, i]
  x = as.matrix(sig[, -1])
  k  = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  para  = rbind (para, k)
  reg   = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  inter = intersect(inter, reg$outlier_detect)
  nout = c(nout, reg$number_outlier)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  beta.fardeep[i, ]   = reg$beta[-1]
  cat('i =', i, '\n')
}

out_gene = NULL
for (i in 1:ncol(outlier)){
  out_tmp = sig$Gene.symbol[outlier[, i]]
  out_gene = c(out_gene, paste0(out_tmp, collapse = ", "))
}
out_dat = data.frame(Samples = colnames(mic), Number.of.outlier = nout, Outlier.genes = out_gene)
write.csv(out_dat, "number_outlier_OV.csv", row.names = FALSE)


library(survival)
library(survminer)
library(RTCGA.clinical)
targ = apply(beta.fardeep[, c(4, 10:12, 15)], 1, sum)
farcl1 = which(targ >= median(targ))
farcl2 = which(targ < median(targ))
far.clust = matrix(0, 514, 1)
far.clust[farcl1] = "high.immune.score"
far.clust[farcl2] = "low.immune.score"
surviv.far = data.frame(status = t(data.matrix(sur[2, ])), day = apply(data.matrix(sur[3:4, ]), 2, 
                function(x) sum(x, na.rm = TRUE)), clust = far.clust)
fit2 = survfit(Surv(day, vital_status) ~ clust, data = surviv.far)
chi2 = survdiff(Surv(day, vital_status) ~ clust, data = surviv.far)$chisq
p2 = pchisq(chi2, length(levels(surviv.far$clust)) - 1, lower.tail = FALSE)
ggsurvplot(fit2, data = surviv.far, pval = TRUE)


################################## CIBERSORT
cib.rela = read.table("Data/CIBERSORT.Output.txt", header = TRUE, row.names = 1, sep = "\t")
cib.score = cib.rela[, 1:22]
targc = apply(cib.score[, c(4, 10:12, 15)], 1, sum)

cibercl1 = which(targc >= median(targc))
cibercl2 = which(targc < median(targc))
ciber.clust = matrix(0, 514, 1)
ciber.clust[cibercl1] = "high.immune.score"
ciber.clust[cibercl2] = "low.immune.score"
surviv.ciber = data.frame(status = t(data.matrix(sur[2, ])), day = apply(data.matrix(sur[3:4, ]), 2, function(x) sum(x, na.rm = TRUE)), clust = ciber.clust)
fit3 = survfit(Surv(day, vital_status) ~ clust, data = surviv.ciber)
chi3 = survdiff(Surv(day, vital_status) ~ clust, data = surviv.ciber)$chisq
p3 = pchisq(chi3, length(levels(surviv.ciber$clust)) - 1, lower.tail = FALSE)
ggsurv3 = ggsurvplot(fit3, data = surviv.ciber, pval = TRUE)

########################## ESTIMATE
library(utils)
rforge = "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(package="estimate")

OvarianCancerExpr = "Data/OV.txt"
filterCommonGenes(input.f = OvarianCancerExpr, output.f="OV_10412genes.gct", id = "GeneSymbol")
estimateScore("OV_10412genes.gct", "OV_estimate_score.gct", platform="affymetrix")
#plotPurity(scores="OV_estimate_score.gct", samples="all_samples", platform="affymetrix")
score = read.table("OV_estimate_score.gct", skip = 2, header = TRUE)
score.tmp1 = score[, -c(1,2)]
score.tmp2 = score.tmp1[, -c(grep("02", temp42), grep("11", temp42))]
score.tmp3 = score.tmp2[, temp6 %in% temp2]
sco = score.tmp3[, order(colnames(score.tmp3))]
rownames(sco) = score[, 1]

scocl1 = which(as.numeric(sco[2,]) > median(as.numeric(sco[2,])))
scocl2 = which(as.numeric(sco[2,]) < median(as.numeric(sco[2,])))
clust = matrix(0, 514, 1)
clust[scocl1] = "high.immune.score"
clust[scocl2] = "low.immune.score"
surviv.est = data.frame(status = t(data.matrix(sur[2, ])), day = apply(data.matrix(sur[3:4, ]), 2, function(x) sum(x, na.rm = TRUE)), clust = clust)
fit1 = survfit(Surv(day, vital_status) ~ clust, data = surviv.est)
chi1 = survdiff(Surv(day, vital_status) ~ clust, data = surviv.est)$chisq
p1 = pchisq(chi1, length(levels(surviv.est$clust)) - 1, lower.tail = FALSE)
ggsurvplot(fit1, data = surviv.est, pval = TRUE)


########### CIBERSORT(not normalize to 1)
source("CIBERSORT/CIBERSORT_no_last_normalization.R")
mix = 2 ^ mic
ciberm = cbind(row.names(mix), mix)
write.table(ciberm, "mixture_file.txt", sep="\t", row.names = FALSE)
write.table(sig, "sig_matrix.txt", sep="\t", row.names = FALSE)

result = CIBERSORT.n(sig_matrix = "./sig_matrix.txt" , 
                     mixture_file = "./mixture_file.txt", perm = 0, QN = FALSE)

cibn.score = result[, 1:22]
id.temp = which(cibn.score < 0)
cibn.score[id.temp] = 0
targcn = apply(cibn.score[, c(4, 10:12, 15)], 1, sum)
cibncl1 = which(targcn > median(targcn))
cibncl2 = which(targcn < median(targcn))
cibn.clust = matrix(0, 514, 1)
cibn.clust[cibncl1] = "high.immune.score"
cibn.clust[cibncl2] = "low.immune.score"
surviv.cibn = data.frame(status = t(data.matrix(sur[2, ])), day = apply(data.matrix(sur[3:4, ]), 2,
                function(x) sum(x, na.rm = TRUE)), clust = cibn.clust)
fit4 = survfit(Surv(day, vital_status) ~ clust, data = surviv.cibn)
chi4 = survdiff(Surv(day, vital_status) ~ clust, data = surviv.cibn)$chisq
p4 = pchisq(chi4, length(levels(surviv.cibn$clust)) - 1, lower.tail = FALSE)
ggsurv4 = ggsurvplot(fit4, data = surviv.cibn, pval = TRUE)



######### whole
whol = rbind(surviv.est, surviv.far, surviv.ciber, surviv.cibn)
Method = c(rep("ESTIMATE", 514), rep("FARDEEP", 514), rep("CIBERSORT", 514), rep("CIBERSORT.modified", 514))
whole = cbind(whol, Method)
whole$Method = factor(whole$Method, levels = c("ESTIMATE", "CIBERSORT", "CIBERSORT.modified", "FARDEEP"))
fit.whole = survfit(Surv(day, vital_status) ~ clust, data = whole)
ggsurvp = ggsurvplot_facet(fit.whole, data = whole, facet.by = "Method", pval = TRUE, pval.size = 10, palette = "jco", censor = FALSE)
ov = ggsurvp + theme_grey() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 15), axis.title=element_text(size=15), axis.text=element_text(size=14), strip.text = element_text(size=15))
save(ov, file = "Figure7_FigureA2/OV/OV.RData")






####################### CIBERSORT abs
Yq = 2 ^ mic
#write.table(Yq[,1:100], "split1.txt", sep="\t", row.names=TRUE, col.names = TRUE)
#write.table(Yq[,101:200], "split2.txt", sep="\t", row.names=TRUE, col.names = TRUE)
#write.table(Yq[,201:300], "split3.txt", sep="\t", row.names=TRUE, col.names = TRUE)
#write.table(Yq[,301:400], "split4.txt", sep="\t", row.names=TRUE, col.names = TRUE)
#write.table(Yq[,401:514], "split5.txt", sep="\t", row.names=TRUE, col.names = TRUE)
cib.abs1 = read.table("Data/absolute1.txt", header = TRUE, row.names = 1, sep = "\t")[, 1:22]
cib.abs2 = read.table("Data/absolute2.txt", header = TRUE, row.names = 1, sep = "\t")[, 1:22]
cib.abs3 = read.table("Data/absolute3.txt", header = TRUE, row.names = 1, sep = "\t")[, 1:22]
cib.abs4 = read.table("Data/absolute4.txt", header = TRUE, row.names = 1, sep = "\t")[, 1:22]
cib.abs5 = read.table("Data/absolute5.txt", header = TRUE, row.names = 1, sep = "\t")[, 1:22]
cib.abs = rbind(cib.abs1,cib.abs2,cib.abs3,cib.abs4,cib.abs5)


targca = apply(cib.abs[, c(4, 10:12, 15)], 1, sum)
cibabscl1 = which(targca > median(targca))
cibabscl2 = which(targca < median(targca))
cibabs.clust = matrix(0, 514, 1)
cibabs.clust[cibabscl1] = "high.immune.score"
cibabs.clust[cibabscl2] = "low.immune.score"
surviv.cibabs = data.frame(status = t(data.matrix(sur[2, ])), day = apply(data.matrix(sur[3:4, ]), 2, function(x) sum(x, na.rm = TRUE)), clust = cibabs.clust)
fit4 = survfit(Surv(day, vital_status) ~ clust, data = surviv.cibabs)
chi4 = survdiff(Surv(day, vital_status) ~ clust, data = surviv.cibabs)$chisq
p4 = pchisq(chi4, length(levels(surviv.cibabs$clust)) - 1, lower.tail = FALSE)
ggsurv4 = ggsurvplot(fit4, data = surviv.cibabs, pval = TRUE, ggtheme = theme_grey())


whol = surviv.cibabs
Method = rep("CIBERSORT.abs", 514)
whole = cbind(whol, Method)
fit.whole = survfit(Surv(day, vital_status) ~ clust, data = whole)
ggsurvp = ggsurvplot_facet(fit.whole, data = whole, facet.by = "Method", pval = TRUE, pval.size = 10, palette = "jco", censor = FALSE)
ov_abs = ggsurvp + theme_grey() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 15), axis.title=element_text(size=15), axis.text=element_text(size=14), strip.text = element_text(size=15))
save(ov_abs, file = "Figure7_FigureA2/OV/OV_abs.RData")

  
  
  
  
  
  