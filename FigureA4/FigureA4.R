lm22 = read.table("Data/LM22.txt", header = TRUE, sep = "\t")
amount = as.matrix(c(2, 1.5, 1.2, 1, 0.8, 0.6, 0.5, 0.3, 0.2))
sig = lm22[, 2:10]
mix = as.matrix(sig) %*% amount
n = nrow(sig)
z = log(sd(as.vector(mix)),2)
mixture = NULL
for (i in 1:50){
  set.seed(i)
  e = 2 ^ rnorm (n, 0, z * 0.3)
  y = mix + e
  mixture = cbind(mixture, y)
}
library(preprocessCore)
mix.quantile = normalize.quantiles(as.matrix(mixture))

sig.c = data.frame(gene = lm22[, 1], cell = sig)
mixture.c = data.frame(gene = lm22[, 1], samp = mix.quantile)
source("CIBERSORT/CIBERSORT_no_last_normalization.R")
source("CIBERSORT/cibersort_without_normalization.R")
write.table(mixture.c, "mixture_file.txt", sep="\t", row.names = FALSE)
write.table(sig.c, "sig_matrix.txt", sep="\t", row.names = FALSE)
m1 = CIBERSORT.nn (sig_matrix = "sig_matrix.txt", mixture_file = "mixture_file.txt", perm=0, QN=F)
m2 = CIBERSORT.n (sig_matrix = "sig_matrix.txt", mixture_file = "mixture_file.txt", perm=0, QN=F)

beta.o = as.vector(m1[, 1:9])
beta.n = as.vector(m2[, 1:9])
beta.t = rep(as.vector(rep(amount, each = 50)), 2)
sse.o = round(sum((beta.o-beta.t)^2), 2)
sse.n = round(sum((beta.n-beta.t)^2), 2)
val1 = c(beta.o, beta.n)
coe = rep(rep(c("B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper", "T.cells.regulatory.Tregs"), each = 50), 2)
meth = rep(c(paste("SVR (SSE = ", sse.o, ")", sep = ""), paste("SVR.znormalization (SSE = ", sse.n, ")", sep = "")), each = 450)
library(ggplot2)


dat = data.frame(Estimated.value = val1, Method = meth, Coefficient = coe)


ref = data.frame(int = amount, slope = 0, Coefficient = c("B.cells.naive", "B.cells.memory", "Plasma.cells", "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.CD4.memory.activated", "T.cells.follicular.helper",  "T.cells.regulatory.Tregs"))

p1 = ggplot(data = dat, aes(x = Method, y = Estimated.value)) + geom_boxplot(aes(fill = Method)) + labs(x = "", y = "Estimated value") + scale_x_discrete(breaks = NULL) +
  facet_wrap( ~ Coefficient, scales="free") + theme_gray() + geom_abline(data = ref, aes(intercept=int, slope = slope), color="red", linetype = "dashed") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold.italic"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 30), axis.title=element_text(size=20), axis.text=element_text(size=20), strip.text = element_text(size=16)) 
