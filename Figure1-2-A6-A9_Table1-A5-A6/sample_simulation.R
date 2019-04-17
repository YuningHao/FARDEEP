sample.sim = function (n = 500, p = 20, a1 = 0.2, a2 = 0.2, noise, outlier, seed){
  set.seed(seed)
  u     = matrix (runif (n * p, 0, 20), n, p)
  rou   = 0.5
  sigma = matrix (rou, p, p)
  diag (sigma) = 1
  x = u %*% (sigma ^ (1/2))
  if (noise == "i"){
    e = rnorm(n, 0, 1)
  }else if (noise == "ii"){
    e = rt(n, 3)
  }
  tau   = matrix(0, n, 1)
  loc.a1  = sample(n, floor(a1 * n))
  if (outlier == "a"){
    tmp = rt(length(loc.a1), 1, ncp = 30)
  }else if (outlier == "b"){
    tmp1 = apply(x[loc.a1, ], 1, max) * 10
    tmp = c()
    for (j in 1:length(tmp1)){
      tmp = c(tmp, rt(1, 1, ncp = tmp1[j]))
      }
    }
  tau[loc.a1] = tmp
  loc.a2  = sample (loc.a1, a2 * length(loc.a1))
  lev = 2 * max(x)
  x[loc.a2, ] = abs(rnorm(length(loc.a2) * p, lev, 1))
  beta  = runif (p, 0, 1)
  y     = x %*% beta + tau + e
  loc   = loc.a1
  result  = list (y = y, x = x, loc = loc, beta = beta)
  return(result)
}



erelated.sim = function(n = 500, p = 20, a1 = 0.2, a2 = 0.2, seed){
  set.seed(seed)
  u     = matrix (runif (n * p, 0, 20), n, p)
  rou   = 0.5
  sigma = matrix (rou, p, p)
  diag (sigma) = 1
  x     = u %*% (sigma ^ (1/2))
  omega = diag(n)
  omega[1:20, 1:20] = 0.7
  omega[21:40, 21:40] = 0.5
  diag(omega) = 1
  e     = matrix(rmvt(1, omega, 3), n, 1)
  tau   = matrix(0, n, 1)
  loc.a1  = sample(n, floor(a1 * n))
  tau[loc.a1] = rt(length(loc.a1), 1, ncp = 30)
  loc.a2  = sample(loc.a1, a2 * length(loc.a1))
  lev = 2 * max(x)
  x[loc.a2, ] = abs(rnorm(length(loc.a2) * p, lev, 1))
  beta  = runif (p, 0, 1)
  y     = x %*% beta + tau + e
  loc   = loc.a1
  result  = list (y = y, x = x, loc = loc, beta = beta)
  return(result)
}
