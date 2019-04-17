source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
library(preprocessCore)
library(affy)
#mixture = NULL
#for (i in 279601 : 279612){
#  cel = ReadAffy(filenames = paste("Data/GSE11058_RAW/GSM", i, ".CEL.gz", sep = ""))
#  e   = exprs(mas5(cel))
#  mixture = cbind(mixture, e)
#}
#cel1 = ReadAffy(filenames = "Data/GSM269529.CEL.gz")
#e1   = exprs(mas5(cel1))
#cel2 = ReadAffy(filenames = "Data/GSM269530.CEL.gz")
#e2   = exprs(mas5(cel2))
#temp1 = data.frame(probe = rownames(e1), tum = (e1 + e2)/2)
#temp2 = data.frame(probe = rownames(mixture), samp = mixture)
#df = merge(temp1, temp2, by = "probe")
#qn = normalize.quantiles(as.matrix(df[, -1]))
#exp = data.frame(probe = df$probe, tum = qn[, 1], mixA = apply(qn[, 2:4], 1, mean),
#                 mixB = apply(qn[, 5:7], 1, mean), mixC = apply(qn[, 8:10], 1, mean),
#                 mixD = apply(qn[, 11:13], 1, mean))
#write.csv(exp, "Figure4/add_tum.csv", row.names = FALSE)

exp = read.csv("Figure4/add_tum.csv")
sig = read.csv("Figure4/GSE11103_matrix_classes.GSE11103-GSE10650.AbbasPure.mas5.bm.K999.0.txt", sep = "")
mix = merge(exp, data.frame(probe = sig$NAME), by = "probe")
z = log(sd(as.matrix(exp[, -1])), 2)

mix_tum = NULL
for (i in c(0, 0.3, 0.6, 0.9)){
  for(j in c(0, 0.3, 0.6, 0.9)){
    tm = j * mix[, 2] + (1 - j) * mix[, 3:6]
    for(w in 1:4){
      set.seed(100 * (i + 1) * j * w)
      up = sample(584, 292)
      down = c(1:584)[-up]
      tm_w = tm[, w]
      tm_w[up] = tm_w[up] + 2 ^ rnorm(292, sd = z * i)
      tm_w[down] = tm_w[down] - 2 ^ rnorm(292, sd = z * i)
      tm_w[tm_w < 0] = 0
      mix_tum = cbind(mix_tum, tm_w)
    }
  }
}

mixA   = c(0.25, 0.125, 0.25, 0.375)
mixB   = c(0.05, 0.317, 0.475, 0.158)
mixC   = c(0.01, 0.495, 0.165, 0.33)
mixD   = c(0.002, 0.333, 0.333, 0.333)


A1 = rbind(mixA, mixB, mixC, mixD)
A2 = rbind(A1, 0.7 * A1, 0.4 *A1, 0.1 * A1)
true.beta = rbind(A2, A2, A2, A2)

source("CIBERSORT/CIBERSORT_no_last_normalization.R")
m = mix_tum
s = sig[, -1]
library(ggplot2)

n.col  = ncol (m)
row_n  = nrow (m)
p      = 4
beta.fardeep = matrix(0, n.col, p)
para  = NULL
nout = NULL
outlier = matrix(0, row_n, n.col)
for (i in 1 : n.col){
  y = as.matrix(m[, i])
  x = as.matrix(s)
  k       = tuningBIC(x = x, y = y, n = row_n, p = p, intercept = TRUE)
  para    = rbind (para, k)
  reg     = fardeep(x = x,  y = y, k = k, intercept = TRUE)
  nout    = c(nout, reg$number_outlier)
  outlier[1 : length(reg$outlier_detect), i] = reg$outlier_detect
  coe     = reg$beta[-1]
  beta.fardeep[i, ] = coe
}
diff.fardeep = abs(beta.fardeep - true.beta)
sse.fardeep  = apply (diff.fardeep, 1, function(x) sum (x^2))

rownames(m) = mix$probe
rownames(s) = sig$NAME 
write.table(m, "mixture_file.txt", sep="\t")
write.table(s, "sig_matrix.txt", sep="\t")
result = CIBERSORT.n (sig_matrix = "sig_matrix.txt" , 
                      mixture_file = "mixture_file.txt", perm=0, QN=F)
beta.ciber = result[ , 1:4]
i = which(beta.ciber < 0)
beta.ciber[i] = 0
diff.ciber = abs(beta.ciber - true.beta)
sse.ciber  = apply (diff.ciber, 1, function(x) sum (x^2))

cat = rep(c("mixA", "mixB", "mixC", "mixD"), 32)
sse.tu = t(cbind(t(sse.fardeep), t(sse.ciber)))
dat.tu = data.frame(Method = factor(rep(c("FARDEEP", "CIBERSORT"), each = 64)), MIX = factor(cat), SSE = sse.tu)
dat.tu$Method = factor(dat.tu$Method, levels=c("CIBERSORT", "FARDEEP"))
dat.tu$MIX= factor(dat.tu$MIX, levels=c("mixA", "mixB", "mixC", "mixD"))
colo = c("darkorchid4", "red")
p.tu1 = ggplot(data = dat.tu, aes(x = MIX, y = SSE)) + geom_boxplot(aes(fill = Method)) + labs(x = "") + scale_fill_manual(values=colo) +
  theme_gray() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 20), axis.title=element_text(size=20), axis.text=element_text(size=20), strip.text = element_text(size=15))

af = as.vector(beta.fardeep)
ac = as.vector(beta.ciber)
at = as.vector(true.beta)
r2.fardeep = round(1 - sum((af-at)^2)/sum((at - mean(at))^2), 2)
r2.ciber = round(1 - sum((ac-at)^2)/sum((at - mean(at))^2), 2)
r.fardeep = round(cor(af, at), 2)
r.ciber = round(cor(ac, at), 2)
coe.tru = rep(at, 2)

coe.est = c(af, ac)
method  = rep(c(paste("FARDEEP (R2 = ",r2.fardeep, ", R = ", r.fardeep, ")", sep = ""), 
                paste("CIBERSORT (R2 = ",r2.ciber, ", R = ", r.ciber, ")", sep = "")), each = 256)
cell.type = rep(rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 64), 2)
dat = data.frame(Method = factor(method), Cell.type = cell.type, coe_est = coe.est, coe_tru = coe.tru)
dat$Method = factor(dat$Method, levels=c(paste("CIBERSORT (R2 = ",r2.ciber, ", R = ", r.ciber, ")", sep = ""), 
                                         paste("FARDEEP (R2 = ",r2.fardeep, ", R = ", r.fardeep, ")", sep = "")))
library(ggplot2)
p.tu2 = ggplot(dat, aes(coe_est, coe_tru, shape = Cell.type, colour = Cell.type)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 3) + labs(y = "True cell amount", x = "Estimated cell amount") + geom_abline(intercept = 0, slope = 1) + 
  facet_wrap( ~ Method, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 20), axis.title=element_text(size=20), axis.text=element_text(size=20), strip.text = element_text(size=16))


### BY CELL

cf = round(c(cor(beta.fardeep[, 1], true.beta[, 1]), cor(beta.fardeep[, 2], true.beta[, 2]), cor(beta.fardeep[, 3], true.beta[, 3]), cor(beta.fardeep[, 4], true.beta[, 4])), 2)
cc = round(c(cor(beta.ciber[, 1], true.beta[, 1]), cor(beta.ciber[, 2], true.beta[, 2]), cor(beta.ciber[, 3], true.beta[, 3]), cor(beta.ciber[, 4], true.beta[, 4])), 2)
sf = round(c(sum((beta.fardeep[, 1] - true.beta[, 1])^2), sum((beta.fardeep[, 2] - true.beta[, 2])^2), sum((beta.fardeep[, 3] - true.beta[, 3])^2), sum((beta.fardeep[, 4] - true.beta[, 4])^2)), 2)
sc = round(c(sum((beta.ciber[, 1] - true.beta[, 1])^2), sum((beta.ciber[, 2] - true.beta[, 2])^2), sum((beta.ciber[, 3] - true.beta[, 3])^2), sum((beta.ciber[, 4] - true.beta[, 4])^2)), 2)
r2f = round(c(1 - sum((beta.fardeep[, 1] - true.beta[, 1])^2)/sum((true.beta[, 1] - mean(true.beta[, 1]))^2), 
              1 - sum((beta.fardeep[, 2] - true.beta[, 2])^2)/sum((true.beta[, 2] - mean(true.beta[, 2]))^2), 
              1 - sum((beta.fardeep[, 3] - true.beta[, 3])^2)/sum((true.beta[, 3] - mean(true.beta[, 3]))^2),
              1 - sum((beta.fardeep[, 4] - true.beta[, 4])^2)/sum((true.beta[, 4] - mean(true.beta[, 4]))^2)), 2)
r2c = round(c(1 - sum((beta.ciber[, 1] - true.beta[, 1])^2)/sum((true.beta[, 1] - mean(true.beta[, 1]))^2), 
              1 - sum((beta.ciber[, 2] - true.beta[, 2])^2)/sum((true.beta[, 2] - mean(true.beta[, 2]))^2), 
              1 - sum((beta.ciber[, 3] - true.beta[, 3])^2)/sum((true.beta[, 3] - mean(true.beta[, 3]))^2),
              1 - sum((beta.ciber[, 4] - true.beta[, 4])^2)/sum((true.beta[, 4] - mean(true.beta[, 4]))^2)), 2)


val = c(af, ac)
tru = rep(at, 2)
meth = rep(c("FARDEEP", "CIBERSORT"), each = 256)
cell = rep(rep(c("Jurkat", "IM-9", "Raji", "THP-1"), each = 64), 2)
dat.r = data.frame(Method = meth, Cell.type = cell, coe_est = val, coe_tru = tru)
library(ggplot2)
library(grid)
p = ggplot(dat.r, aes(coe_est, coe_tru, colour = Method)) + xlim(0, 1) + ylim(0, 1) +
  geom_point(cex = 2) + labs(y = "True cell amount", x = "Estimated cell fraction") + geom_abline(intercept = 0, slope = 1) + scale_colour_manual(values=colo) +
  facet_wrap( ~ Cell.type, scales="free") + theme_gray() + theme(plot.title = element_text(hjust = 0.8), legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 18), axis.title=element_text(size=18), strip.text = element_text(size=18))



text1 = paste("FARDEEP \n (R = ", cf[1], ", SSE = ", sf[1], ", R2 = ", r2f[1], ") \n  CIBERSORT \n (R = ", cc[1], ", SSE = ", sc[1], ", R2 = ", r2c[1], ")", sep = "")
text2 = paste("FARDEEP \n (R = ", cf[2], ", SSE = ", sf[2], ", R2 = ", r2f[2], ")  \n CIBERSORT \n (R = ", cc[2], ", SSE = ", sc[2], ", R2 = ", r2c[2], ")", sep = "")
text3 = paste("FARDEEP \n (R = ", cf[3], ", SSE = ", sf[3], ", R2 = ", r2f[3], ")  \n CIBERSORT \n (R = ", cc[3], ", SSE = ", sc[3], ", R2 = ", r2c[3], ")", sep = "")
text4 = paste("FARDEEP \n (R = ", cf[4], ", SSE = ", sf[4], ", R2 = ", r2f[4], ")  \n CIBERSORT \n (R = ", cc[4], ", SSE = ", sc[4], ", R2 = ", r2c[4], ")", sep = "")

text = c(text1, text2, text3, text4)
cell = c("Jurkat", "IM-9", "Raji", "THP-1")
graphLabels = data.frame(Cell.type = cell, text = text)


bycell = p + geom_text(data = graphLabels, aes(x = 0.25, y = 0.9, label = text, color = NULL), size = 3)



library(cowplot)
library(PASWR)
p.tu3 = plot_grid(p.tu1, p.tu2, bycell, ncol = 1, labels = c("A", "B", "C"), label_size = 20)


