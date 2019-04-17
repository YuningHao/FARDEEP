sample.sim = function (n = 500, p = 20, sig = 1, a1 = 0.2, a2 = 0.2){
  ### n : sample size
  ### p : number of predictors
  ### sig : variance of y
  ### a1  : proportion of outliers
  ### a2  : proportion of leverage points in outliers
  ### lev : mean of leverage ponits on x 
  u     = matrix (runif (n * p, 0, 20), n, p)
  rou   = 0.5
  sigma = matrix (rou, p, p)
  diag (sigma) = 1
  x     = u %*% (sigma ^ (1/2))
  e     = rnorm (n, 0, sig)
  tau   = matrix(0, n, 1)
  loc.a1  = sample (n, floor(a1 * n))
  tau[loc.a1] = rnorm (length(loc.a1), 30, 5)
  loc.a2  = sample (loc.a1, a2 * length(loc.a1))
  lev = 2 * max(x)
  x[loc.a2, ] = rnorm(length(loc.a2) * p, lev, 1)
  beta  = runif (p, 0, 1)
  y     = x %*% beta + tau + e
  loc   = loc.a1
  result  = list (y = y, x = x, loc = loc, beta = beta)
  return(result)
}

source("sourcecode/fardeep_function.R")
source("sourcecode/Tuning_BIC.R")
set.seed(100)
sam = sample.sim (n = 500, p = 20, sig = 1, a1 = 0.2, a2 = 0.2)
x   = sam$x
y   = sam$y
beta = sam$beta
loc  = sam$loc
write.table (x, file= "fardeepx.txt", row.names=F, col.names=F, quote=F)
write.table (y, file= "fardeepy.txt", row.names=F, col.names=F, quote=F)
write.table (beta, file= "fardeepbeta.txt", row.names=F, col.names=F, quote=F)
write.table (loc, file= "fardeeploc.txt", row.names=F, col.names=F, quote=F)


n = nrow(x)
p = ncol(x)
O = 0.2 * 500
library(nnls)

tp.fardeep1 = NULL
fp.fardeep1 = NULL
number.outlier1 = NULL
for (i in seq (0.1, 0.5, 0.05)) {
  tun1   = tuningBIC(x = x, y = y, n = n, p = p, alpha1 = i, alpha2 = 1.5, up = 10, low = 1 , nn = TRUE, intercept = TRUE, lognorm = FALSE)
  reg1  = fardeep (x = x,  y = y, alpha1 = i, alpha2 = 1.5, k = tun1, nn = TRUE, intercept = TRUE)  
  out1   = reg1$outlier_detect
  numo1  = reg1$number_outlier
  tpr1   = 1 - sum (is.na(match (loc, out1))) / O
  fpr1   = sum (is.na(match (out1, loc))) / (n - O)
  tp.fardeep1 = c (tp.fardeep1, tpr1)
  fp.fardeep1 = c (fp.fardeep1, fpr1)
  number.outlier1 = c (number.outlier1, numo1)
}
tp.fardeep2 = NULL
fp.fardeep2 = NULL
number.outlier2 = NULL
for (i in seq (1.1, 1.9, 0.1)) {
  tun2   = tuningBIC(x = x,  y = y, n = n, p = p, alpha1 = 0.1, alpha2 = i, up = 10, low = 1 , nn = TRUE, intercept = FALSE, lognorm = FALSE)
  reg2   = fardeep (x = x,  y = y, alpha1 = 0.1, alpha2 = i, k = tun2, nn = TRUE, intercept = FALSE)  
  out2   = reg2$outlier_detect
  numo2  = reg2$number_outlier
  tpr2   = 1 - sum (is.na(match (loc, out2))) / O
  fpr2   = sum (is.na(match (out2, loc))) / (n - O)
  tp.fardeep2 = c (tp.fardeep2, tpr2)
  fp.fardeep2 = c (fp.fardeep2, fpr2)
  number.outlier2 = c (number.outlier2, numo2)
}

tp.fardeep4 = NULL
fp.fardeep4 = NULL
number.outlier4 = NULL
for (i in seq (1, 10, 0.1)) {
  reg4   = fardeep (x = x,  y = y, alpha1 = 0.1, alpha2 = 1.5, k = i, nn = TRUE, intercept = TRUE)  
  tun4   = tuningBIC(x = x,  y = y, n = n, p = p, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1 , nn = TRUE, intercept = TRUE, lognorm = FALSE)
  out4   = reg4$outlier_detect
  numo4  = reg4$number_outlier
  tpr4   = 1 - sum (is.na(match (loc, out4))) / O
  fpr4   = sum (is.na(match (out4, loc))) / (n - O)
  tp.fardeep4 = c (tp.fardeep4, tpr4)
  fp.fardeep4 = c (fp.fardeep4, fpr4)
  number.outlier4 = c (number.outlier4, numo4)
}
