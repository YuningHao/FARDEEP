# library(annotate)
# library(AnnotationDbi)
# library(affy)
# sig = read.table("Data/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt",
#                  header = TRUE) 
# cel1 = ReadAffy(filenames = "Data/GSE32474_RAW/GSM803626_113395hp133a11.cel.gz")
# e1   = exprs(mas5(cel1))
# cel2 = ReadAffy(filenames = "Data/GSE32474_RAW/GSM803685_113454hp133a11.cel.gz")
# e2   = exprs(mas5(cel2))
# cel3 = ReadAffy(filenames = "Data/GSE32474_RAW/GSM803744_118180hp133a11.cel.gz")
# e3   = exprs(mas5(cel3))
# ee2 = as.data.frame((e1 + e2 + e3)/2)
# ee2$NAME = rownames(ee2)
# dat = merge(sig, ee2, by = "NAME")
# write.csv(ee2, "Data/MDA_MB_435.csv", row.names = FALSE)

MDA_MB_435 = read.csv("Data/MDA_MB_435.csv")
sig = read.table("Data/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt",
                 header = TRUE) 
dat = merge(sig, MDA_MB_435, by = "NAME")
mix1 = as.matrix(dat[, -1]) %*% matrix(rep(0.2, 5), 5, 1)
mix2 = as.matrix(dat[, -1]) %*% matrix(c(rep(0.05, 4), 0.8), 5, 1)
mix = data.frame(mix1 = mix1, mix2 = mix2)
rownames(mix) = dat$NAME

## CIBERSORT
source("CIBERSORT/CIBERSORT.R")
s = sig[, -c(1, 6)]
row.names(s) = dat$NAME
write.table(mix, "mixture_file.txt", sep="\t")
write.table(s, "sig_matrix.txt", sep="\t")
result = CIBERSORT(sig_matrix = "sig_matrix.txt" , 
                   mixture_file = "mixture_file.txt", perm = 0, QN = F)
beta.c = result[ , 1:4]

## FARDEEP
source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
n = nrow(mix)
p = 4
n.col = 2
beta.f = matrix(0, n.col, p)
outlier = matrix(0, n, n.col)
for (i in 1 : n.col){
  y = mix[, i]
  x = as.matrix(sig[, -1])
  k = tuningBIC(x = x, y = y, n = n, p = p, intercept = TRUE)
  reg = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  beta.f[i, ]   = reg$beta[-1]
}

beta.f.relative = beta.f/apply(beta.f, 1, sum)  # FARDEEP relative


## outlier analysis
MDA_MB_435 = dat$GSM803626_113395hp133a11.cel.gz
order_id = order(MDA_MB_435)

outlier1 = outlier[, 1]
outlier_id1 = which(order_id %in% outlier1)
plot(MDA_MB_435[order_id], pch = 20, xlab = "Index", ylab = "Expression value", axes = FALSE, main = "MDA-MB-435: MIXTURE 1")
points(outlier_id1, MDA_MB_435[order_id[outlier_id1]], col = "red", pch = 20)
axis(1)
axis(2)

outlier2 = outlier[, 2]
outlier_id2 = which(order_id %in% outlier2)
plot(MDA_MB_435[order_id], pch = 20, xlab = "Index", ylab = "Expression value", axes = FALSE, main = "MDA-MB-435: MIXTURE 2")
points(outlier_id2, MDA_MB_435[order_id[outlier_id2]], col = "red", pch = 20)
axis(1)
axis(2)
