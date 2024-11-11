library(monotone)
Phi_Durot <- function(Y,alpha,sd,c){
  n = length(Y)
  F_n = cumsum(Y)
  S_n = -Inf
  x = seq(1/n,1,1/n)
  for (l in 1:(n-1)) {
    for (i in 1:(n-l)) {
      data = Y[i:(i+l)]
      est = legacy(data,number = 2)
      gcm = cumsum(est)
      gcm = gcm - max(gcm-F_n[i:(i+l)])
      S_n = max(max(sqrt(1/((l+1)*n))/sd*(F_n[i:(i+l)]-gcm)),S_n)
      if(S_n>=c){
        return(1)
      }
    }
  }
  
  return(0)
}

Simulate_c = function(n,alpha,R){
  S_n = rep(-Inf,R)
  for (r in 1:R) {
    Y = rnorm(n,0,1)
    F_n = cumsum(Y)
    for (l in 1:(n-1)) {
      for (i in 1:(n-l)) {
        data = Y[i:(i+l)]
        est = legacy(data,number = 2)
        gcm = cumsum(est)
        gcm = gcm - max(gcm-F_n[i:(i+l)])
        S_n[r] = max(max(sqrt(1/((l+1)*n))*(F_n[i:(i+l)]-gcm)),S_n[r])
      }
    }
  }
  return(quantile(S_n,1-alpha))
}




