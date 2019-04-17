# install.packages("hgu133plus2hsentrezgcdf_18.0.0.tar.gz",type = "source",repos=NULL)
# download.file("http://mbni.org/customcdf/18.0.0/entrezg.download/hgu133plus2hsentrezg.db_18.0.0.tar.gz", method = "auto",destfile = "hgu133plus2hsentrezg.db_18.0.0.tar.gz")
# install.packages("hgu133plus2hsentrezg.db_18.0.0.tar.gz",type = "source",repos=NULL)
# download.file(url = "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/18.0.0/affy_1.40.0.tar.gz", method = "auto",destfile = "affy_1.40.0.tar.gz")
# install.packages("affy_1.40.0.tar.gz",type = "source",repos=NULL)
# library(affy)
# library(hgu133plus2hsentrezgcdf)
# library(hgu133plus2hsentrezg.db)
# library(preprocessCore)
# cel = ReadAffy(cdfname="hgu133plus2hsentrezgcdf")
# e   = exprs(mas5(cel))
# fl = normalize.quantiles(as.matrix(e))
# ex = data.frame(entrez = gsub(row.names(e), pattern="_at",replacement=""), ex = fl)
# cols = c("entrez", colnames(e))
# colnames(ex) = cols
# columns(hgu133plus2hsentrezg.db)
# OUT = as.data.frame(mapIds(hgu133plus2hsentrezg.db, keys = as.character(ex$entrez), column = "SYMBOL", keytype = "ENTREZID"))
# colnames(OUT) = "gene"
# data = data.frame(entrez = ex$entrez, gene = OUT$gene)
# ge = na.omit(merge(data, ex, by = "entrez"))[, -1]
# unique(ge$gene)
# write.csv(ge, "GSE65135_Follicular_Lymphoma.csv", row.names = FALSE)

setwd("Figure5/Follicular_lymphoma")
ge = read.csv("GSE65135_Follicular_Lymphoma.csv")
lm22 = read.table("../../Data/LM22.txt", header = TRUE, sep = "\t")
name.s = lm22$Gene.symbol
name.m = ge$gene
gep = ge[name.m %in% name.s, ]
idm = order(gep$gene)
gexp = gep[idm, -1]
row.names(gexp) = gep$gene[idm]
sg = lm22[name.s %in% name.m, ]
ids = order(sg$Gene.symbol)
sig = sg[ids, ]

####################fardeep
source("../../sourcecode/fardeep_function.R")
source("../../sourcecode/Tuning_BIC.R")
source("../../CIBERSORT/CIBERSORT.R")
n     = nrow(gexp)
p     = ncol(sig) - 1
n.col = ncol(gexp)
beta.fardeep = matrix(0, n.col, p)
betaabs.fardeep = matrix(0, n.col, p)
para  = NULL
outlier = matrix(0, n, n.col)
nout = NULL
inter = 1:n
for (i in 1 : n.col){
  y = gexp[, i]
  x = as.matrix(sig[, -1])
  k = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  para  = rbind (para, k)
  reg   = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  inter = intersect(inter, reg$outlier_detect)
  nout  = c(nout, reg$number_outlier)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  beta.fardeep[i, ] = reg$beta[-1]/sum(reg$beta[-1])
  betaabs.fardeep[i, ]  = reg$beta[-1]
}

inter_gene = lm22$Gene.symbol[inter]
out_gene = NULL
for (i in 1:ncol(outlier)){
  out_tmp = sig$Gene.symbol[outlier[, i]]
  out_gene = c(out_gene, paste0(out_tmp, collapse = ", "))
}
out_dat = data.frame(Samples = paste("GSM", 1587831 : 1587844, sep = ""), Number.of.outlier = nout, Outlier.genes = out_gene)
write.csv(out_dat, "number_outlier_fl.csv", row.names = FALSE)

#"CD79A"  "BCL2A1"

######### CIBERSORT
ematrix = ge
write.table(ematrix, file="FL.txt", sep="\t", col.names=T, row.names=F, quote=FALSE)
ciber.result = read.table("./CIBERSORT.Output_relative.txt", header = TRUE, sep = "\t")
beta.ciber = ciber.result[, 2:23]


# normalization
b.f = apply(beta.fardeep[, 1:2], 1, sum)
b.c = apply(beta.ciber[, 1:2], 1, sum)
cd4.f = apply(beta.fardeep[, 5:9], 1, sum)
cd4.c = apply(beta.ciber[, 5:9], 1, sum)
cd8.f = beta.fardeep[, 4]
cd8.c = beta.ciber[, 4]

b.fardeep = b.f/(b.f + cd4.f + cd8.f)
b.ciber = b.c/(b.c + cd4.c + cd8.c)
cd4.fardeep = cd4.f/(b.f + cd4.f + cd8.f)
cd4.ciber = cd4.c/(b.c + cd4.c + cd8.c)
cd8.fardeep = cd8.f/(b.f + cd4.f + cd8.f)
cd8.ciber = cd8.c/(b.c + cd4.c + cd8.c)

############ true
cd4.true = c(0.1735, 0.1939, 0.3333, 0.2247, 0.0882, 0.4023, 0.2143, 0.1368, 0.0753, 0.1122, 0.1262, 0.0891, 0.2759, 0.0505)
cd8.true = c(0.0204, 0.0612, 0.069, 0.0225, 0.0196, 0.0575, 0.0714, 0.0316, 0.0108, 0.0102, 0.0291, 0.0297, 0.0575, 0.0202)
b.true = c(0.8061, 0.7449, 0.5977, 0.7528, 0.8922, 0.5402, 0.7143, 0.8316, 0.9104, 0.8776, 0.8447, 0.8812, 0.6667, 0.9293)

true = c(cd4.true, cd8.true, b.true)
ciber = c(cd4.ciber, cd8.ciber, b.ciber)
fard = c(cd4.fardeep, cd8.fardeep, b.fardeep)

cor.f = round(cor(fard, true, method = "pearson"), 2)
cor.c = round(cor(ciber, true, method = "pearson"), 2)
sse.f = round(sum((fard - true)^2), 2)
sse.c = round(sum((ciber - true)^2), 2)
r2.f = round(1 - sum((fard - true)^2)/sum((true - mean(true))^2), 3)
r2.c = round(1 - sum((ciber - true)^2)/sum((true - mean(true))^2), 3)

library(ggplot2)
coe.tru = rep(true, 2)
coe.est = c(fard, ciber)
method  = rep(c(paste("FARDEEP (R2 = ", r2.f, ", R = ", cor.f, ", SSE = ", sse.f, ")", sep = ""), paste("CIBERSORT (R2 = ", r2.c, ", R = ", cor.c, ", SSE = ", sse.c, ")", sep = "")), each = 42)
cell.type = rep(rep(c("CD4 T cells", "CD8 T cells", "B cells"), each = 14), 2)
dat = data.frame(Method = factor(method), Cell.type = cell.type, coe_est = coe.est, coe_tru = coe.tru)
dat$Method = factor(dat$Method, levels = c(paste("CIBERSORT (R2 = ", r2.c, ", R = ", cor.c, ", SSE = ", sse.c, ")", sep = ""), paste("FARDEEP (R2 = ", r2.f, ", R = ", cor.f, ", SSE = ", sse.f, ")", sep = "")))
flfig = ggplot(dat, aes(coe_est, coe_tru, shape = Cell.type, colour = Cell.type)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(y = "Flow cytometry fraction", x = "Estimated cell fraction") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Method, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))

#ggsave("FL.tiff", plot = flfig, width = 40, height = 30, units = "cm", dpi = 350)

save(flfig, file = "flfig.RData")



################################## by cell
cf = round(c(cor(b.fardeep, b.true), cor(cd4.fardeep, cd4.true), cor(cd8.fardeep, cd8.true)), 2)
cc = round(c(cor(b.ciber, b.true), cor(cd4.ciber, cd4.true), cor(cd8.ciber, cd8.true)), 2)
sf = round(c(sum((b.fardeep - b.true)^2), sum((cd4.fardeep - cd4.true)^2), sum((cd8.fardeep - cd8.true)^2)), 2)
sc = round(c(sum((b.ciber - b.true)^2), sum((cd4.ciber - cd4.true)^2), sum((cd8.ciber - cd8.true)^2)), 2)
r2f = round(c(1 - sum((b.fardeep - b.true)^2)/sum((b.true - mean(b.true))^2), 1 - sum((cd4.fardeep - cd4.true)^2)/sum((cd4.true - mean(cd4.true))^2), 1 - sum((cd8.fardeep - cd8.true)^2)/sum((cd8.true - mean(cd8.true))^2)), 2)
r2c = round(c(1 - sum((b.ciber - b.true)^2)/sum((b.true - mean(b.true))^2), 1 - sum((cd4.ciber - cd4.true)^2)/sum((cd4.true - mean(cd4.true))^2), 1 - sum((cd8.ciber - cd8.true)^2)/sum((cd8.true - mean(cd8.true))^2)), 2)


val = c(b.fardeep, cd4.fardeep, cd8.fardeep, b.ciber, cd4.ciber, cd8.ciber)
tru = c(b.true, cd4.true, cd8.true)
meth = rep(c("FARDEEP", "CIBERSORT"), each = 42)
cell = rep(rep(c("B cells", "CD4 T cells", "CD8 T cells"), each = 14), 2)
dat.r = data.frame(Method = meth, Cell.type = cell, coe_est = val, coe_tru = tru)
library(ggplot2)
library(grid)
p = ggplot(dat.r, aes(coe_est, coe_tru, colour = Method)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(y = "Flow cytometry fraction", x = "Estimated cell fraction") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Cell.type, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))



text1 = paste("FARDEEP \n (R = ", cf[1], ", SSE = ", sf[1], ", R2 = ", r2f[1], ") \n      CIBERSORT \n (R = ", cc[1], ", SSE = ", sc[1], ", R2 = ", r2c[1], ")", sep = "")
text2 = paste(" FARDEEP \n (R = ", cf[2], ", SSE = ", sf[2], ", R2 = ", r2f[2], ")  \n CIBERSORT \n (R = ", cc[2], ", SSE = ", sc[2], ", R2 = ", r2c[2], ")", sep = "")
text3 = paste("FARDEEP \n (R = ", cf[3], ", SSE = ", sf[3], ", R2 = ", r2f[3], ")  \n CIBERSORT \n (R = ", cc[3], ", SSE = ", sc[3], ", R2 = ", r2c[3], ")", sep = "")
text = c(text1, text2, text3)
cell = c("B cells", "CD4 T cells", "CD8 T cells")
graphLabels = data.frame(Cell.type = cell, text = text)

by.cell = p + geom_text(data = graphLabels, aes(x = 0.4, y = 0.9, label = text, color = NULL), size = 5)

save(by.cell, file = "fl_bycell.RData")
