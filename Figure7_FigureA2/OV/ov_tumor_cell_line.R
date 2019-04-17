ov = read.csv("~/Dropbox/Yuning/PLOS/number_outlier_ov.csv")
lusc = read.csv("~/Dropbox/Yuning/PLOS/number_outlier_lusc.csv")

ov.gene = unlist(sapply(ov$Outlier.genes, function(x) unlist(strsplit(as.character(x), ", "))))
freq.ov = as.data.frame(ftable(ov.gene))
freq.ov = freq.ov[order(freq.ov$Freq, decreasing = TRUE),]

write.csv(freq.ov, "~/Dropbox/Yuning/PLOS/freq_ov.csv")

lusc.gene = unlist(sapply(lusc$Outlier.genes, function(x) unlist(strsplit(as.character(x), ", "))))
freq.lusc = as.data.frame(ftable(lusc.gene))
freq.lusc = freq.lusc[order(freq.lusc$Freq, decreasing = TRUE),]

write.csv(freq.lusc, "~/Dropbox/Yuning/PLOS/freq_lusc.csv")


library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)
library(affy)
library(preprocessCore)
setwd("Data/OV_cell_line/")
cel = ReadAffy(cdfname="hgu133plus2hsentrezgcdf")
exp = exprs(mas5(cel))
ex = data.frame(entrez = gsub(row.names(exp), pattern="_at",replacement=""), ex = normalize.quantiles(as.matrix(exp)))
cols = c("entrez", colnames(exp))
colnames(ex) = cols
OUT = as.data.frame(mapIds(hgu133plus2hsentrezg.db, keys = as.character(ex$entrez), column = "SYMBOL", keytype = "ENTREZID"))
colnames(OUT) = "gene"
data = data.frame(entrez = ex$entrez, gene = OUT$gene)
ge = na.omit(merge(data, ex, by = "entrez"))[, -1]
lm22 = read.table("../LM22.txt", header = TRUE, sep = "\t")
df = merge(data.frame(gene = lm22$Gene.symbol), ge, by = "gene")
NCI_ADR_RES = data.frame(gene = df$gene, NCI_ADR_RES = round(apply(df[, c(2, 9, 16)], 1, mean)))
OVCAR_3 = data.frame(gene = df$gene, OVCAR_3 = round(apply(df[, c(4, 11, 18)], 1, mean)))
OVCAR_4 = data.frame(gene = df$gene, OVCAR_4 = round(apply(df[, c(5, 12, 19)], 1, mean)))
OVCAR_5 = data.frame(gene = df$gene, OVCAR_5 = round(apply(df[, c(6, 13, 20)], 1, mean)))
OVCAR_8 = data.frame(gene = df$gene, OVCAR_8 = round(apply(df[, c(7, 14, 21)], 1, mean)))
SK_OV_3 = data.frame(gene = df$gene, SK_OV_3 = round(apply(df[, c(8, 15, 22)], 1, mean)))
IGROV1 = data.frame(gene = df$gene, IGROV1 = round(apply(df[, c(3, 10, 17)], 1, mean)))

ov.cell.line = freq.ov[1:22, ]
colnames(ov.cell.line) = c("gene", "removal Frequency (n = 514)")
ov.cell.line = merge(ov.cell.line, NCI_ADR_RES, by = "gene")
ov.cell.line = merge(ov.cell.line, OVCAR_3, by = "gene")
ov.cell.line = merge(ov.cell.line, OVCAR_4, by = "gene")
ov.cell.line = merge(ov.cell.line, OVCAR_5, by = "gene")
ov.cell.line = merge(ov.cell.line, OVCAR_8, by = "gene")
ov.cell.line = merge(ov.cell.line, SK_OV_3, by = "gene")
ov.cell.line = merge(ov.cell.line, IGROV1, by = "gene")
ov.cell.line = ov.cell.line[order(ov.cell.line$`removal Frequency (n = 514)`, decreasing = TRUE),] 
write.csv(ov.cell.line, "~/Dropbox/Yuning/PLOS/Figure6_FigureS2/OV/ov_cell_line.csv", row.names = FALSE)

