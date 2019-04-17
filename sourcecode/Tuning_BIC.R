tuningBIC = function (x, y, n, p, alpha1 = 0.1, alpha2 = 1.5, up = 10, low = 1, nn = TRUE, 
                      intercept = FALSE, lognorm = TRUE, step.size = 0.1){
  ##################################################################################################
  ########################## Tuning Parameter using adjusted BIC ###################################
  ##################################################################################################
  ### INPUT
  ### x : independent variables
  ### y : dependent variables
  ### n : sample size
  ### p : number of independent variables
  ### up : upper bound for k
  ### low : lower bound for k
  ### alpha1 : between 0 and 1
  ### alpha2 : larger than 1
  ### nn : non-negative coefficients
  ### intercept : include intercept in the model
  ### lognorm : whether noise is log-normal distributed, default TRUE
  ###################################
  ### OUTPUT
  ### k : parameter k for FARDEEP
  para  = NULL
  BIC.fardeep = NULL
  for (j in seq (low, up, step.size)) {
    reg1  = fardeep(x = x, y = y, alpha1 = alpha1, alpha2 = alpha2, k = j, nn = nn, intercept = intercept)
    if (intercept){
      if (lognorm){
        res = log(abs(reg1$Y.new - cbind (1, reg1$X.new) %*% reg1$beta), 2)
      }else{
        res = reg1$Y.new - cbind (1, reg1$X.new) %*% reg1$beta
      }
    }else{
      if (lognorm){
        res = log(abs(reg1$Y.new - reg1$X.new %*% reg1$beta), 2)
      }else{
        res = reg1$Y.new - reg1$X.new %*% reg1$beta
      }
    }
    sse   = t (res) %*% res
    t     = reg1$number_outlier + p + 1
    no    = reg1$number_outlier
    BIC2  = (n - no) * log (sse / (n - no)) + t * (log(n - no) + 1)
    BIC.fardeep = rbind (BIC.fardeep, BIC2) 
    cat("k =", j, "  BIC:", BIC.fardeep)
  }
  ind2    = which.min (BIC.fardeep)
  seq_dat     = seq (low, up, step.size)     
  k       = seq_dat[ind2]
  return (k)
}