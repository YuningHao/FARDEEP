sample.sim = function (n, p, sig, a1, a2, nneg = TRUE, intercept = FALSE){
  ### n : sample size
  ### p : number of predictors
  ### sig : variance of y
  ### a1  : proportion of outliers
  ### a2  : proportion of leverage points in outliers
  ### nneg : non-negtive or not
  ### intercept : including intercept or not
  u     = matrix (runif (n * p, 0, 30), n, p)
  sigma = matrix (0.5, p, p)
  diag (sigma) = 1
  x     = u %*% (sigma ^ (1/2))
  e     = rnorm (n, 0, sig)
  tau   = matrix(0, n, 1)
  loc.a1  = sample (n, floor(a1 * n))
  tau[loc.a1] = rnorm (length(loc.a1), 20, 5)
  loc.a2  = sample (loc.a1, a2 * length(loc.a1))
  lev = 2 * max(x)
  x[loc.a2, ] = rnorm(length(loc.a2) * p, lev, 1)
  if (nneg){
    if (intercept){
      beta  = runif(p + 1, 0, 1)
      y     = cbind(1, x) %*% beta + tau + e
    }else{
      beta  = runif(p, 0, 1)
      y     = x %*% beta + tau + e
      }
  }else{
    if(intercept){
      beta  = rnorm(p + 1, 0, 1)
      y     = cbind(1, x) %*% beta + tau + e
    }else{
      beta  = rnorm(p, 0, 1)
      y     = x %*% beta + tau + e
    }
  }
  loc   = loc.a1
  result  = list (y = y, x = x, loc = loc, beta = beta)
  return(result)
}

##### Sample with 20 % outliers including 20% leverage point (non-negative)
##### generate sample with outlier
sam1 = sample.sim (n = 200, p = 5, sig = 1, a1 = 0.2, a2 = 0.2)
##### tuning parameter k2
k = tuningBIC (sam1$x, sam1$y, 200, 5, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1)
##### ALTS
m1 = alts(sam1$x,  sam1$y, alpha1 = 0.2, alpha2 = 1.1, k = k)


##### Sample with 10 % outliers including 5% leverage point
##### generate sample with outlier
sam2 = sample.sim (n = 300, p = 10, sig = 1, a1 = 0.1, a2 = 0.05, nneg = FALSE)
##### tuning parameter k2
k = tuningBIC (sam2$x, sam2$y, 300, 10, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = FALSE, 
               intercept = FALSE)
##### ALTS
m2 = alts(sam2$x,  sam2$y, alpha1 = 0.2, alpha2 = 1.1, k = k, nn = FALSE, intercept = FALSE)



##### Sample with 15 % outliers including 10% leverage point (with intercept)
##### generate sample with outlier
sam3 = sample.sim (n = 500, p = 15, sig = 1, a1 = 0.15, a2 = 0.1, nneg = FALSE, intercept = TRUE)
##### tuning parameter k2
k = tuningBIC (sam3$x, sam3$y, 500, 15, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = FALSE, 
               intercept = TRUE)
##### ALTS
m3 = alts(sam3$x,  sam3$y, alpha1 = 0.1, alpha2 = 1.5, k = k, nn = FALSE, intercept = TRUE)

