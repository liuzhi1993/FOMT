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

support = function(x,n,h){
  res_1 = max(floor_FOMT((x-h)*n),1)
  res_2 = min(ceiling_FOMT((x+h)*n),n)
  return(c(res_1,res_2))
}

Q = function(x_1,x_2,x,h,kernel = "epanechnikov", k=0){
  n1 = length(x_1)
  n2 = length(x_2)
  res = matrix(1,n1,n2)
  if(k==1){
    for (i in 1:n1) {
      for (j in 1:n2) {
        res[i,j] = abs(x_1[i]-x_2[j])
      }
    }
  }
  temp1 = K((x_1-x)/h,kernel)
  temp2 = K((x_2-x)/h,kernel)
  temp = temp1%*%t(temp2)
  res = res*temp
  return(res)
}

K = function(x,kernel){
  if(kernel == "epanechnikov"){
    res = 3/4*(1-x^2)*(abs(x)<=1)
  }
  if(kernel == "triangle"){
    res = (1-abs(x))*(abs(x)<=1)
  }
  if(kernel == "quartic"){
    res = (15/16*(1-x^2)^2)*(abs(x)<=1)
  }
  if(kernel == "triweight"){
    res = (35/32*(1-x^2)^3)*(abs(x)<=1)
  }
  if(kernel == "cosinus"){
    res = (pi/4*cos(pi*x/2))*(abs(x)<=1)
  }
  return(res)
}

test_C = function(X,Y,x,h,sd = 1,kernel = "epanechnikov",k=0){
  n = length(Y)
  supp = support(x,n,h)
  i1 = supp[1]
  i2 = supp[2]
  Y = Y[(X<=i2/n)&(X>=i1/n)]
  X = X[(X<=i2/n)&(X>=i1/n)]
  temp = Q(X,X,x,h,kernel,k)
  diff_Y = matrix(0,length(Y),length(Y))
  sign_X = matrix(0,length(X),length(X))
  for (i in 1:length(Y)) {
    for (j in 1:length(Y)) {
      diff_Y[i,j] = Y[i]-Y[j]
      sign_X[i,j] = sign(X[j]-X[i])
    }
  }
  res1 = 1/2*diff_Y*sign_X*temp
  res1 = sum(res1)
  
  temp2 = sign_X*temp
  temp2 = colSums(temp2) 
  temp2 = temp2*temp2
  res2 = sum(temp2)
  res2 = sd*sqrt(res2)
  res = res1/res2
  return(res)
}

simulate_c <- function(n,X = NA,R = 100, kernel="epanechnikov",k=0, alpha = 0.05){
  u = 0.5
  if((length(X)==1) && is.na(X)){
    X = seq(1/n,1,1/n)
    h_max = (n-1)/(2*n)
  }else{
    X = runif(n,0,1)
    h_max = 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        h_max = max(abs(X[i]-X[j]),h_max) 
      }
    }
    h_max = h_max/2
  }
  h_min = 2/5*h_max*(log(n)/n)^(1/3)
  l = floor(log2(5/2*(n/log(n))^(1/3)))
  l = seq(0,l,1)
  H_n = h_max*u^l
  c = rep(0,R)
  for (r in 1:R) {
    #print(r)
    error = rnorm(n,0,1)
    test_statistic = test_C(X,error,X[1],h_max,1,kernel,k)
    for(i in 1:n){
      for (j in 1:length(H_n)) {
        test_statistic = max(test_C(X,error,X[i],H_n[j],1,kernel,k),test_statistic)
      }
    }
    c[r] = test_statistic
  }
  critical_value = quantile(c,1-alpha)
  print(critical_value)
  return(critical_value)
}

Phi_C <- function(X,Y,kernel="epanechnikov",sd =1,k=0,c){
  n = length(Y)
  u = 0.5
  h_max = 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      h_max = max(abs(X[i]-X[j]),h_max) 
    }
  }
  h_max = h_max/2
  h_min = 2/5*h_max*(log(n)/n)^(1/3)
  l = floor(log2(5/2*(n/log(n))^(1/3)))
  l = seq(0,l,1)
  H_n = h_max*u^l
  critical_value = c
  test_statistic = test_C(X,Y,X[1],h_max,sd,kernel,k)
  for(i in 1:n){
    for (j in 1:length(H_n)) {
      test_statistic = max(test_C(X,Y,X[i],H_n[j],sd,kernel,k),test_statistic)
      if(test_statistic >= critical_value){
        #print(1)
        return(1)
      }
    }
  }
  #print(test_statistic)
  if(test_statistic >= critical_value){
    #print(1)
    return(1)
  }else{
    #print(0)
    return(0)
  }
}

