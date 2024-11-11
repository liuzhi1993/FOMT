support = function(x,n,h){
  res_1 = max(floor_FOMT((x-h)*n),1)
  res_2 = min(ceiling_FOMT((x+h)*n),n)
  return(c(res_1,res_2))
}

C = function(r){
  return(sqrt(2*log(1/r)))
}

ceiling_FOMT = function(x){
  res = ceiling(x)
  if(x-res==0){
    res = res + 1
  }
  return(res)
}

floor_FOMT = function(x){
  res = floor(x)
  if(x-res==0){
    res = res -1
  }
  return(res)
}

psi = function(x, type = 1){
  if(type ==1 ){
    res = (abs(x)<=1)*(1-abs(x))
  }else{
    res = (abs(x)<=1)*(1-abs(x))*x
    #res = (abs(x)<=2)*(1-abs(x)/2)*x
  }
  return(res)
}

PPsi = function(Y,t,h,sd,type=1){
  n = length(Y)
  supp = support(t,n,h)
  Y = Y[supp[1]:supp[2]]
  X = seq(supp[1],supp[2],1)/n
  X = (X-t)/h
  temp = psi(X,type)
  if(sum(temp*temp)==0){
    return(0)
  }else{
    res = sum(temp*Y)/(sqrt(sum(temp*temp))*sd)
  }
  return(res)
}

Simulate_c = function(n, R = 100, type = 1){
  S_n = seq(1/n,1/2,1/n)
  sd = 1
  test_statistic = rep(-Inf,R)
  for (r in 1:R) {
    data = rnorm(n,0,sd)
    for (i in 1:length(S_n)) {
      h = S_n[i]
      L_n = seq(h,1-h,h)
      temp1 = rep(0,length(L_n))
      for(j in 1:length(L_n)){
        t = L_n[j]
        temp1[j] = PPsi(data,t,h,sd,type)
      }
      #print(temp1)
      if(type ==1){
        if(length(L_n)>=2){
          for (k in 1:((length(L_n)-1))) {
            for (l in (k+1):length(L_n)) {
              test_statistic[r] = max(test_statistic[r], (temp1[l]-temp1[k])/2-C(2*h))
            }
          }
        }
      }else{
        test_statistic[r] = max(test_statistic[r],max(temp1)-C(2*h))
      }
    }
    #print(c(r,test_statistic[r]))
  }
  kappa = quantile(test_statistic,1-alpha)
}

Phi_DS = function(Y,sd,type=1,critical_value){
  n = length(Y)
  S_n = seq(1/n,1/2,1/n)
  test_statistic = -Inf
  for (i in 1:length(S_n)) {
    h = S_n[i]
    L_n = seq(h,1-h,h)
    temp1 = rep(0,length(L_n))
    for(j in 1:length(L_n)){
      t = L_n[j]
      temp1[j] = PPsi(Y,t,h,sd,type)
    }
    if(type ==1){
      if(length(L_n)>=2){
        for (k in 1:((length(L_n)-1))) {
          for (l in (k+1):length(L_n)) {
            test_statistic = max(test_statistic, (temp1[l]-temp1[k])/2-C(2*h))
            if(test_statistic>=critical_value){
              return(1)
            }
          }
        }
      }
    }else{
      test_statistic = max(test_statistic,max(temp1)-C(2*h))
      if(test_statistic>=critical_value){
        return(1)
      }
    }
  }
  return(0)
}

