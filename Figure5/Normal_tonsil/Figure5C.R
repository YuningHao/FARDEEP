# setwd("Figure5/Normal_tonsil")
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
# ttq = normalize.quantiles(as.matrix(e))
# ex = data.frame(entrez = gsub(row.names(e), pattern="_at",replacement=""), ex = ttq)
# cols = c("entrez", colnames(e))
# colnames(ex) = cols
# columns(hgu133plus2hsentrezg.db)
# OUT = as.data.frame(mapIds(hgu133plus2hsentrezg.db, keys = as.character(ex$entrez), column = "SYMBOL", keytype = "ENTREZID"))
# colnames(OUT) = "gene"
# data = data.frame(entrez = ex$entrez, gene = OUT$gene)
# ge = na.omit(merge(data, ex, by = "entrez"))[, -1]
# unique(ge$gene)
# write.csv(ge, "GSE65135_Normal_Tonsil.csv", row.names = FALSE)


setwd("Figure5/Normal_tonsil")
ge = read.csv("GSE65135_Normal_Tonsil.csv")
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
para  = NULL
nout = NULL
outlier = matrix(0, n, n.col)
for (i in 1 : n.col){
  y = gexp[, i]
  x = as.matrix(sig[, -1])
  k  = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  para  = rbind (para, k)
  reg   = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  nout = c(nout, reg$number_outlier)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  beta.fardeep[i, ]   = reg$beta[-1]/sum(reg$beta[-1])
}

out_gene = NULL
for (i in 1:ncol(outlier)){
  out_tmp = sig$Gene.symbol[outlier[, i]]
  out_gene = c(out_gene, paste0(out_tmp, collapse = ", "))
}
out_dat = data.frame(Samples = paste("GSM", 1587845 : 1587854, sep = ""), Number.of.outlier = nout, Outlier.genes = out_gene)
write.csv(out_dat, "number_outlier_tonsil.csv", row.names = FALSE)

######### CIBERSORT
ematrix = ge
write.table(ematrix, file="tonsil.txt", sep="\t", col.names=T, row.names=F, quote=FALSE)
ciber.result = read.table("CIBERSORT.Output_relative.txt", header = TRUE, sep = "\t")
beta.ciber = ciber.result[, 2:23]



b.f = apply(beta.fardeep[, 1:2], 1, sum)
b.c = apply(beta.ciber[, 1:2], 1, sum)
t.f = apply(beta.fardeep[, 4:10], 1, sum)
t.c = apply(beta.ciber[, 4:10], 1, sum)
o.f = apply(beta.fardeep[, c(3, 11:22)], 1, sum)
o.c = apply(beta.ciber[, c(3, 11:22)], 1, sum)



library(ggplot2)

cell = rep(rep(c("B cells", "T cells", "Other"), each = 10), 2)
meth = rep(c("FARDEEP", "CIBERSORT"), each = 30)
val = c(b.f, t.f, o.f, b.c, t.c, o.c)
samp = rep(c(paste("B", c(1, 3:6), sep = ""), paste("T", c(1, 3:6), sep = "")), 6)
bar = data.frame(cell = cell, meth = meth, Proportion = val, Sample = samp)
bar$cell = factor(bar$cell, levels= c("Other", "T cells", "B cells"))
bar$meth = factor(bar$meth, levels= c("CIBERSORT", "FARDEEP"))

tonsil = ggplot(data = bar, aes(x = Sample, y = Proportion, fill = cell)) +
  geom_bar(stat = "identity") + facet_wrap(~ meth, nrow = 2) + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title.x = element_blank(), axis.title.y =element_text(size=18), axis.text = element_text(size=15), strip.text = element_text(size=18))

ggsave("tonsil.tiff", plot = tonsil, width = 40, height = 30, units = "cm", dpi = 350)
save(tonsil, file = "tonsil.RData")

