# setwd("Figure5/GSE20300/GSE20300_RAW/")
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
# ex = data.frame(entrez = gsub(row.names(e), pattern="_at",replacement=""), ex = normalize.quantiles(as.matrix(e)))
# cols = c("entrez", colnames(e))
# colnames(ex) = cols
# columns(hgu133plus2hsentrezg.db)
# OUT = as.data.frame(mapIds(hgu133plus2hsentrezg.db, keys = as.character(ex$entrez), column = "SYMBOL", keytype = "ENTREZID"))
# colnames(OUT) = "gene"
# data = data.frame(entrez = ex$entrez, gene = OUT$gene)
# ge = na.omit(merge(data, ex, by = "entrez"))[, -1]
# unique(ge$gene)
# write.csv(ge, "../GSE20300.csv", row.names = FALSE)


setwd("Figure5/GSE20300")
ge = read.csv("GSE20300.csv")
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

out_gene = NULL
for (i in 1:ncol(outlier)){
  out_tmp = sig$Gene.symbol[outlier[, i]]
  out_gene = c(out_gene, paste0(out_tmp, collapse = ", "))
}
out_dat = data.frame(Samples = paste("GSM", 508883 : 508906, sep = ""), Number.of.outlier = nout, Outlier.genes = out_gene)
write.csv(out_dat, "number_outlier_GSE20300.csv", row.names = FALSE)



source("../../CIBERSORT/CIBERSORT.R")
s = sig[, -1]
m = gexp
row.names(s) = sig$Gene.symbol
write.table(m, "mixture_file.txt", sep="\t")
write.table(s, "sig_matrix.txt", sep="\t")
result = CIBERSORT(sig_matrix = "sig_matrix.txt" , 
                   mixture_file = "mixture_file.txt", perm = 0, QN = F)
beta.ciber = result[ , 1:22]

neu = c(76.4, 38, 70.7, 62.1, 90.5, 65.8, 39.2, 58.5, 60.5, 89.4, 43.6, 91.3, 55.4, 14.5, 81, 63.4, 50.1, 54.8, 47.7, 53.8, 59.9, 77, 30.4, 47.8)
lym = c(19.4, 47.5, 17.1, 20.3, 2, 25, 37.1, 28, 27.8, 6.3, 44, 1.6, 31.3, 61.5, 17.9, 23.1, 30.6, 15.2, 37.8, 37.2, 32.6, 18.5, 52, 39.7)
mono = c(3.6, 12.4, 11.8, 14.6, 4.8, 7.3, 14.1, 9.6, 8.3, 3.7, 5.8, 4.6, 8.4, 21.9, 0.4, 11.6, 18.1, 27.8, 12.1, 8.4, 5.8, 4.1, 9.4, 6.6)
eos = c(0.4, 1.4, 0.1, 2.5, 2.4, 1.4, 9.2, 3.3, 2.8, 0.1, 6.5, 2.5, 4.6, 1.3, 0.6, 1.6, 1.2, 1.6, 1.6, 0.4, 1.5, 0.2, 7.7, 4.7)
baso = c(0.2, 0.8, 0.3, 0.5, 0.3, 0.5, 0.4, 0.6, 0.6, 0.5, 0, 0, 0.3, 0.8, 0.1, 0.3, 0, 0.6, 0.8, 0.2, 0.3, 0.2, 0.5, 1.2)

truth = data.frame(Neutrophils = neu, Lymphocytes = lym, Monocytes = mono, Eosinophils = eos)/100

beta.f = data.frame(Neutrophils = beta.fardeep[, 22], Lymphocytes = apply(beta.fardeep[, 1:12], 1, sum),
                    Monocytes = beta.fardeep[, 13], Eosinophils = beta.fardeep[, 21])
beta.c = data.frame(Neutrophils = beta.ciber[, 22], Lymphocytes = apply(beta.ciber[, 1:12], 1, sum),
                    Monocytes = beta.ciber[, 13], Eosinophils = beta.ciber[, 21])

tr.prop = truth/apply(truth, 1, sum)

f.prop = beta.f/apply(beta.f, 1, sum)

c.prop = beta.c/apply(beta.c, 1, sum)


f.all = as.vector(as.matrix(f.prop))
c.all = as.vector(as.matrix(c.prop))
t.all = as.vector(as.matrix(tr.prop))

cor.f = round(cor(f.all, t.all), 2)
cor.c = round(cor(c.all, t.all), 2)
sse.f = round(sum((f.all - t.all)^2), 2)
sse.c = round(sum((c.all - t.all)^2), 2)
r2.f = round(1 - sum((f.all - t.all)^2)/sum((t.all - mean(t.all))^2), 2)
r2.c = round(1 - sum((c.all - t.all)^2)/sum((t.all - mean(t.all))^2), 2)

library(ggplot2)
coe.tru = rep(t.all, 2)
coe.est = c(f.all, c.all)
method  = rep(c(paste("FARDEEP (R2 = ", r2.f, ", R = ", cor.f, ", SSE = ", sse.f, ")", sep = ""), paste("CIBERSORT (R2 = ", r2.c, ", R = ", cor.c, ", SSE = ", sse.c, ")", sep = "")), each = 96)
cell.type = rep(rep(c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils"), each = 24), 2)
dat = data.frame(Method = factor(method), Cell.type = cell.type, coe_est = coe.est, coe_tru = coe.tru)
dat$Method = factor(dat$Method, levels = c(paste("CIBERSORT (R2 = ", r2.c, ", R = ", cor.c, ", SSE = ", sse.c, ")", sep = ""), paste("FARDEEP (R2 = ", r2.f, ", R = ", cor.f, ", SSE = ", sse.f, ")", sep = "")))
GSE20300 = ggplot(dat, aes(coe_est, coe_tru, shape = Cell.type, colour = Cell.type)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(y = "Coulter counter fractions", x = "Estimated cell fraction") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Method, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))

save(GSE20300, file = "GSE20300.RData")


################################## by cell
cf = round(c(cor(f.prop[, 1], tr.prop[, 1]), cor(f.prop[, 2], tr.prop[, 2]), cor(f.prop[, 3], tr.prop[, 3]), cor(f.prop[, 4], tr.prop[, 4])), 2)
cc = round(c(cor(c.prop[,1], tr.prop[, 1]), cor(c.prop[, 2], tr.prop[, 2]), cor(c.prop[, 3], tr.prop[, 3]), cor(c.prop[, 4], tr.prop[, 4])), 2)
sf = round(c(sum((f.prop[, 1] - tr.prop[, 1])^2), sum((f.prop[, 2] - tr.prop[, 2])^2), sum((f.prop[, 3] - tr.prop[, 3])^2), sum((f.prop[, 4] - tr.prop[, 4])^2)), 2)
sc = round(c(sum((c.prop[, 1] - tr.prop[, 1])^2), sum((c.prop[, 2] - tr.prop[, 2])^2), sum((c.prop[, 3] - tr.prop[, 3])^2), sum((c.prop[, 4] - tr.prop[, 4])^2)), 2)
r2f = round(c(1 - sum((f.prop[, 1] - tr.prop[, 1])^2)/sum((tr.prop[, 1] - mean(tr.prop[, 1]))^2), 1 - sum((f.prop[, 2] - tr.prop[, 2])^2)/sum((tr.prop[, 2] - mean(tr.prop[, 2]))^2), 
              1 - sum((f.prop[, 3] - tr.prop[, 3])^2)/sum((tr.prop[, 3] - mean(tr.prop[, 3]))^2), 1 - sum((f.prop[, 4] - tr.prop[, 4])^2)/sum((tr.prop[, 4] - mean(tr.prop[, 4]))^2)), 2)
r2c = round(c(1 - sum((c.prop[, 1] - tr.prop[, 1])^2)/sum((tr.prop[, 1] - mean(tr.prop[, 1]))^2), 1 - sum((c.prop[, 2] - tr.prop[, 2])^2)/sum((tr.prop[, 2] - mean(tr.prop[, 2]))^2),
              1 - sum((c.prop[, 3] - tr.prop[, 3])^2)/sum((tr.prop[, 3] - mean(tr.prop[, 3]))^2), 1 - sum((c.prop[, 4] - tr.prop[, 4])^2)/sum((tr.prop[, 4] - mean(tr.prop[, 4]))^2)), 2)


val = c(f.all, c.all)
tru = rep(t.all, 2)
meth = rep(c("FARDEEP", "CIBERSORT"), each = 96)
cell = rep(rep(c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils"), each = 24), 2)
dat.r = data.frame(Method = meth, Cell.type = cell, coe_est = val, coe_tru = tru)
library(ggplot2)
library(grid)
p = ggplot(dat.r, aes(coe_est, coe_tru, colour = Method)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(x = "Coulter counter fractions", y = "Estimated cell fraction") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Cell.type, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))



text1 = paste("FARDEEP \n (R = ", cf[1], ", SSE = ", sf[1], ", R2 = ", r2f[1], ") \n  CIBERSORT \n (R = ", cc[1], ", SSE = ", sc[1], ", R2 = ", r2c[1], ")", sep = "")
text2 = paste("FARDEEP \n (R = ", cf[2], ", SSE = ", sf[2], ", R2 = ", r2f[2], ")  \n CIBERSORT \n (R = ", cc[2], ", SSE = ", sc[2], ", R2 = ", r2c[2], ")", sep = "")
text3 = paste("FARDEEP \n (R = ", cf[3], ", SSE = ", sf[3], ", R2 = ", r2f[3], ")  \n CIBERSORT \n (R = ", cc[3], ", SSE = ", sc[3], ", R2 = ", r2c[3], ")", sep = "")
text4 = paste("FARDEEP \n (R = ", cf[4], ", SSE = ", sf[4], ", R2 = ", r2f[4], ")  \n CIBERSORT \n (R = ", cc[4], ", SSE = ", sc[4], ", R2 = ", r2c[4], ")", sep = "")
text = c(text1, text2, text3, text4)
cell = c("Neutrophils", "Lymphocytes", "Monocytes", "Eosinophils")
graphLabels = data.frame(Cell.type = cell, text = text)


gse20300.by.cell = p + geom_text(data = graphLabels, aes(x = 0.4, y = 0.9, label = text, color = NULL), size = 5)

save(gse20300.by.cell, file = "GSE20300_bycell.RData")
