# library(annotate)
# library(AnnotationDbi)
# library(affy)
# sig = read.table("Data/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt",
#                  header = TRUE) 
# cel1 = ReadAffy(filenames = "Data/GSM269529.CEL.gz")
# e1   = exprs(mas5(cel1))
# cel2 = ReadAffy(filenames = "Data/GSM269530.CEL.gz")
# e2   = exprs(mas5(cel2))
# ee2  = as.data.frame((e1 + e2)/2)
# ee2$NAME = rownames(ee2)
# dat = merge(sig, ee2, by = "NAME")
# write.csv(ee2, "Data/HCT116.csv", row.names = FALSE)


HCT116 = read.csv("Data/HCT116.csv")
sig = read.table("Data/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt",
                 header = TRUE) 
dat = merge(sig, HCT116, by = "NAME")
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
hc116 = dat$GSM269529.CEL.gz
order_id = order(hc116)
par(mfrow=c(2, 2))

outlier1 = outlier[, 1]
outlier_id1 = which(order_id %in% outlier1)
plot(hc116[order_id], pch = 20, xlab = "Index", ylab = "Expression value", axes = FALSE, main = "HCT116: MIXTURE 1")
points(outlier_id1, hc116[order_id[outlier_id1]], col = "red", pch = 20)
axis(1)
axis(2)

outlier2 = outlier[, 2]
outlier_id2 = which(order_id %in% outlier2)
plot(hc116[order_id], pch = 20, xlab = "Index", ylab = "Expression value", axes = FALSE, main = "HCT116: MIXTURE 2")
points(outlier_id2, hc116[order_id[outlier_id2]], col = "red", pch = 20)
axis(1)
axis(2)