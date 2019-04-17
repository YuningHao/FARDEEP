fardeep = function (x,  y, alpha1 = 0.1, alpha2 = 1.5, k = 10, nn = TRUE, intercept = FALSE){
  ##################################################################################################
  ###################### FAST And Robust DEconvolution of Expression Profiles ######################
  ##################################################################################################
  ### INPUT
  ### x: Independent variables
  ### y: Dependent variables
  ### alpha1 : between 0 and 1
  ### alpha2 : larger than 1
  ### k : larger than 1
  ### nn : non-negative coefficients
  ### intercept : include intercept in the model
  ##################################################################################################
  ### OUTPUT 
  ### beta : estimation of coefficients
  ### number_outlier : number of outliers
  ### outlier_detect : index of detected outliers
  ### X.new : good observed points for independent variables
  ### Y.new : good observed points for dependent variables
  ### k : modified k (if the input value is not appropriate)
  if(nn){
    library (nnls)
    if(intercept){
      m1 = nnls (cbind(1, x), y)
    }else{
      m1 = nnls (x, y)
    }
  }else{
    if(intercept){
      m1 = lm (y ~ x)
    }else{
      m1 = lm (y ~ x - 1) 
    }
  }
  res = abs (resid(m1))
  n   = length (y)
  order_id = order (res, decreasing = F)
  res_int  = res[order_id]
  Y_int = y[order_id]
  X_int = x[order_id, ]
  index1 = min (which (res_int > median (res_int)))
  k_up   = n - index1
  y_out_up = Y_int[c ((n - k_up + 1) : n)]
  k.low.ex = alpha1 * k_up
  k_low    = ceiling (k.low.ex)
  out_id   = (1:n > (n - k_low))
  kep_id   = (1:n < (n - k_low + 1))
  Y_new    = Y_int[kep_id]
  X_new    = X_int[kep_id, ]
  ii = 1
  
  repeat{
    if(nn){
      if(intercept){
        m_new = nnls (cbind(1, X_new), Y_new)
        beta_new = m_new$x
        Y_hat    = cbind(1, x) %*% beta_new
      }else{
        m_new = nnls (X_new, Y_new)
        beta_new = m_new$x
        Y_hat    = x %*% beta_new
      }
    }else{
      if(intercept){
        m_new = lm (Y_new ~ X_new)
        beta_new = m_new$coefficients
        Y_hat    =  cbind (1, x) %*% beta_new
      }else{
        m_new = lm (Y_new ~ X_new - 1) 
        beta_new = m_new$coefficients
        Y_hat    = x %*% beta_new
      }
    }
    res_new  = abs (y - Y_hat)
    order_id = order (res_new, decreasing = F)
    res_new_ord = res_new[order_id]
    Y_ord    = y[order_id]
    X_ord    = x[order_id, ]
    index1   = min (which (res_new_ord > k * median (res_new_ord)))
    temp     = n - index1
    if(temp <= 0) {
      while (temp <= 0) {
        k = k - 1
        index1 = min (which (res_int > k * median (res_int)))
        temp   = n - index1
      }
    }
    k_up     = min (temp, k_up) # make sure k_up is non-increasing
    k.low.ex = alpha2 * k.low.ex   # update N_low
    k_low    = min (ceiling(k.low.ex), k_up) # prevent the lower bound exceeding the upper bound
    out_id   = (1:n > (n - k_low))
    Y_new    = Y_ord[!out_id]
    X_new    = X_ord[!out_id, ]
    if (k_low >= k_up){
      break
    }
    cat(ii, "\n")
    ii = ii + 1
  }
  if(nn){
    if(intercept){
      model = nnls (cbind(1, X_new), Y_new)
    }else{
      model = nnls (X_new, Y_new)
    }
    coefficients = model$x
  }else{
    if(intercept){
      model = lm (Y_new ~ X_new)
    }else{
      model = lm (Y_new ~ X_new - 1)
    }
    coefficients = model$coefficients
  }
  number_outlier = sum(out_id)
  outlier_id     = order_id[out_id]
  result         = list(beta = coefficients, number_outlier = number_outlier, outlier_detect = outlier_id, 
                        X.new = X_new, Y.new = Y_new, k = k)
  return(result)
}
