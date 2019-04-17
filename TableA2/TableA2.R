lmsig  = read.table("Data/LM22.txt", header = TRUE, sep = "\t")
frac   = NULL
for (i in 1 : 9){
  set.seed(i)
  f = round(runif(22, 0, 1), 2)
  frac = rbind(frac, f)
}
names  = as.matrix (lmsig[ , 1])
lmm    = as.matrix (lmsig[, -1]) %*% t(frac)
lmsig  = lmsig[, -1]
row.names(lmm) = names
row.names(lmsig) = names

source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
source("CIBERSORT/CIBERSORT_no_last_normalization.R")
library(glmnet)
library(nnls)
lm22     = function(percent = 0.1){
  y      = lmm
  z      = log(sd(as.vector(y)), 2)
  x      = as.matrix (lmsig)
  n.col  = ncol (y)
  row_n  = nrow (y)
  p      = ncol(x)
  sampled_n = floor (percent * row_n)
  ind    = sample (row_n, sampled_n)
  ind.down  = ind[1 : floor(sampled_n/2)]
  ind.up    = ind[(floor(sampled_n/2) + 1) : sampled_n]
  down   = sample(row_n, floor(row_n/2))
  num.s  = 1 : row_n
  up     = num.s[-down]
  for (i in 1 : n.col){
    y[down, i]  = y[down, i] - 2 ^ rnorm (length (y[down, i]), mean = 0, sd = z * 0.1)
    y[up, i] = y[up, i] + 2 ^ rnorm (length (y[up, i]), mean = 0, sd = z * 0.1)
    y[ind.down, i]  = y[ind.down, i] - 2 ^ rnorm (length (y[ind.down, i]), mean = 10, sd = z * 0.3)
    y[ind.up, i] = y[ind.up, i] + 2 ^ rnorm (length (y[ind.up, i]), mean = 10, sd = z * 0.3)
    neg = which (y[ , i] < 0)
    y[neg, i] = 0
  }
  write.table(y, file= paste("TableA2/y", percent, ".txt", sep = ""), row.names=F, col.names=F, quote=F)
  write.table(cbind(row.names(y), y), file= paste("TableA2/yciber", percent, ".txt", sep = ""), row.names=F, sep = "\t")
  coe.fardeep = NULL
  coe.nnls = NULL
  coe.dcq = NULL
  R2.fardeep = NULL
  R2.nnls = NULL
  R2.dcq = NULL
  R.fardeep = NULL
  R.nnls = NULL
  R.dcq = NULL
  for (i in 1 : n.col){
    ### fardeep
    y.fardeep   = as.numeric(y[ , i])
    k = tuningBIC(x = x, y = y.fardeep, n = row_n, p = p, intercept = TRUE, lognorm = TRUE)
    reg = fardeep(x = x, y = y.fardeep, k = k, intercept = TRUE)
    af  = reg$beta[-1]
    coe.fardeep = rbind(coe.fardeep, af)
    at = frac[i, ]
    R2temp.fardeep = round(1 - sum((af - at)^2)/sum((at - mean(at))^2), 3)
    R2.fardeep = cbind(R2.fardeep, R2temp.fardeep)
    Rtemp.fardeep = round(cor(af, at), 3)
    R.fardeep = cbind(R.fardeep, Rtemp.fardeep)
    ### nnls
    coe.nn = nnls(x, y.fardeep)$x
    coe.nnls = rbind(coe.nnls, coe.nn)
    R2temp.nnls = round(1 - sum((coe.nn - at)^2)/sum((at - mean(at))^2), 3)
    R2.nnls = cbind(R2.nnls, R2temp.nnls)
    Rtemp.nnls = round(cor(coe.nn, at), 3)
    R.nnls = cbind(R.nnls, Rtemp.nnls)
    ### dcq
    elastic = glmnet (x, y.fardeep, alpha = 0.05, lambda.min.ratio = 0.2, family = c('gaussian'), nlambda = 100,
                      intercept = FALSE)$beta
    coe.d = elastic[, ncol(elastic)]
    coe.dcq = rbind (coe.dcq, coe.d)
    R2temp.dcq = round(1 - sum((coe.d - at)^2)/sum((at - mean(at))^2), 3)
    R2.dcq = cbind(R2.dcq, R2temp.dcq)
    Rtemp.dcq = round(cor(coe.d, at), 3)
    R.dcq = cbind(R.dcq, Rtemp.dcq)
  }
  sse.fardeep = apply(coe.fardeep - frac, 1, function(x) sum (x^2))
  sse.nnls = apply (coe.nnls - frac, 1, function(x) sum (x^2))
  sse.dcq = apply (coe.dcq - frac, 1, function(x) sum (x^2))
  write.table(y, "TableA2/mixture_file.txt", sep="\t")
  write.table(x, "TableA2/sig_matrix.txt", sep="\t")
  beta.cc = CIBERSORT.n(sig_matrix = "TableA2/sig_matrix.txt" ,
                        mixture_file = "TableA2/mixture_file.txt", perm=0, QN=F)[ , 1 : 22]
  index = which (beta.cc < 0)
  beta.cc[index] = 0
  sse.cc = apply (beta.cc - frac, 1, function(x) sum (x^2))
  R2.ciber = NULL
  R.ciber = NULL
  for (j in 1 : 9){
    at = frac[j, ]
    ac = beta.cc[j, ]
    R2temp.ciber = round(1 - sum((ac - at)^2)/sum((at - mean(at))^2), 3)
    R2.ciber = cbind(R2.ciber, R2temp.ciber)
    Rtemp.ciber = round(cor(ac, at), 3)
    R.ciber = cbind(R.ciber, Rtemp.ciber)
  }
  result  = list(out = ind, R2.ciber = R2.ciber, R2.fardeep = R2.fardeep, R2.nnls = R2.nnls, R2.dcq = R2.dcq,
                 sse.ciber = sse.cc, sse.fardeep = sse.fardeep, sse.nnls = sse.nnls, sse.dcq = sse.dcq,
                 R.ciber = R.ciber, R.fardeep = R.fardeep, R.nnls = R.nnls, R.dcq = R.dcq)
  return(result)
}

sse.c = matrix (0, 25, 9)
sse.f = matrix (0, 25, 9)
sse.n = matrix (0, 25, 9)
sse.d = matrix (0, 25, 9)

r2.c = matrix (0, 25, 9)
r2.f = matrix (0, 25, 9)
r2.n = matrix (0, 25, 9)
r2.d = matrix (0, 25, 9)

r.c = matrix (0, 25, 9)
r.f = matrix (0, 25, 9)
r.n = matrix (0, 25, 9)
r.d = matrix (0, 25, 9)

for (k in 1 : 25){
  set.seed(k + 100)
  result = lm22(k/50)
  sse.c[k, ] = result$sse.ciber
  sse.f[k, ] = result$sse.fardeep
  sse.n[k, ] = result$sse.nnls
  sse.d[k, ] = result$sse.dcq
  r2.c[k, ] = result$R2.ciber
  r2.f[k, ] = result$R2.fardeep
  r2.n[k, ] = result$R2.nnls
  r2.d[k, ] = result$R2.dcq
  r.c[k, ] = result$R.ciber
  r.f[k, ] = result$R.fardeep
  r.n[k, ] = result$R.nnls
  r.d[k, ] = result$R.dcq
}

mix = NULL
for (i in 1 : 25){
  y.tmp = as.matrix(read.table(paste("TableA2/y", i/50, ".txt", sep = "")))
  mix = cbind(mix, y.tmp)
}
write.table(mix, "TableA2/mix_for_pert.txt", col.names = FALSE, row.names = FALSE)


coe.pert = read.table("pert/coe_pert_tableA2.txt")
y.true = NULL
for(i in 1 : 25){
  y.true = rbind(y.true, frac)
}

sse.p = t(matrix(apply(coe.pert - y.true, 1, function(x) sum(x ^ 2)), 9, 25))
r.p = c()
for (i in 1:225){
  r.p = c(r.p, round(cor(as.numeric(coe.pert[i, ]), y.true[i,]), 3))
}
r2.p = c()
for (i in 1:225){
  r2.p = c(r2.p, round(1 - sum((as.numeric(coe.pert[i, ]) - y.true[i,]) ^ 2)/sum((y.true[i,] - mean(y.true[i,])) ^ 2), 3))
}

range(sse.c)
range(sse.f)
range(sse.n)
range(sse.d)
range(sse.p)

range(r2.c)
range(r2.f)
range(r2.n)
range(r2.d)
range(r2.p)

range(r.c)
range(r.f)
range(r.n)
range(r.d)
range(r.p)


mean.c = as.matrix(apply(sse.c, 1, mean))
mean.f = as.matrix(apply(sse.f, 1, mean))
mean.n = as.matrix(apply(sse.n, 1, mean))
mean.p = as.matrix(apply(sse.p, 1, mean))
mean.d = as.matrix(apply(sse.d, 1, mean))

apply(sse.c, 1, sd)
as.matrix(apply(sse.f, 1, sd))
apply(sse.n, 1, sd)
apply(sse.p, 1, sd)
apply(sse.d, 1, sd)


