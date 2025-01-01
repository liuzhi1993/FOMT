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

KL_max = function(kernel){
  if(kernel == "epanechnikov"){
    res = c(3/4,3/2)
  }
  if(kernel == "triangle"){
    res = c(1,1)
  }
  if(kernel == "quartic"){
    res = c(15/16,5*sqrt(3)/6)
  }
  if(kernel == "cosinus"){
    res = c(pi/4,pi*pi/8)
  }
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

weight_vector = function(x,n,h,beta,kernel){
  supp = support(x,n,h)
  res = rep(0,supp[2]-supp[1]+1)
  X = seq(supp[1],supp[2],1)/n
  X = (X-x)/h
  if(beta == 1 || ((beta == 2) && (x>=h) && (x<=1-h))){
    temp = K(X,kernel)
    res = temp/sum(temp)
  }else{
    temp = K(X,kernel)
    alpha_0 = sum(temp)
    alpha_1 = sum(X*temp)
    alpha_2 = sum(X*X*temp)
    res = (alpha_2-alpha_1*X)*temp/(alpha_0*alpha_2-alpha_1^2)
  } 
  return(res)
}

LP = function(I,Y,h,beta,kernel,middel_weight){
  n = length(Y)
  supp = support(I/n,n,h)
  y = Y[supp[1]:supp[2]]
  if((h<=I/n)&&(I/n<=1-h)){
    w = middel_weight
  }else{
    w = weight_vector(I/n,n,h,beta,kernel)
  }
  if(length(y)<length(w)){
    d = length(w)-length(y)
    if(d%%2==0){
      q = d/2
      y = c(rep(0,q),y,rep(0,q))
    }else{
      q = (d-1)/2
      y = c(rep(0,q+1),y,rep(0,q))
    }
  }
  if(length(y)>length(w)){
    d = length(y)-length(w)
    if(d%%2==0){
      q = d/2
      w = c(rep(0,q),w,rep(0,q))
    }else{
      q = (d-1)/2
      w = c(rep(0,q+1),w,rep(0,q))
    }
  }
  res = sum(y*w)
  return(res)
}

CALM = function(Y,A,kernel,kappa = 1.5,sd = 0.1, scale_A = 0.5){
  #print(scale_A)
  n = length(Y)
  M = ceiling_FOMT(4/5*log(n,4)+1/5*log(log(n),4)+1)
  est = matrix(NA,M,length(A))
  lambda_0 = 0.5
  if(kernel == "epanechnikov"){
    C_rho = 4*sqrt(3/5)*sd/lambda_0
  }
  if(kernel == "triangle"){
    C_rho = 4*sqrt(2/3)*sd/lambda_0
  }
  if(kernel == "quartic"){
    C_rho = 4*sqrt(5/7)*sd/lambda_0
  }
  if(kernel == "cosinus"){
    C_rho = pi*sd/lambda_0
  }
  C_rho = C_rho * scale_A
  m = 1 
  h = 1/n*4^(m-1)
  middel_weight = weight_vector(0.5, n, h, 1, kernel)
  for (i in 1:length(A)) {
    est[m,i] = LP(A[i],Y,h,1,kernel,middel_weight)
    if(is.na(est[m,i])){
      est[m,i] = Y[i]
    }
  }
  for (m in 2:M) {
    m_bar = m 
    h = 1/n*4^(m-1)
    middel_weight = weight_vector(0.5, n, h, 1, kernel)
    for (i in 1:length(A)) {
      est[m,i] = LP(A[i],Y,h,1,kernel,middel_weight)
    }
    for (k in 1:(m-1)) {
      #print(c(k,m))
      #print(max(abs(est[k,]-est[m,])))
      #print(4*kappa*C_rho*sqrt(log(n)/4^(k-1)))
      if(max(abs(est[k,]-est[m,]))>=4*kappa*C_rho*sqrt(log(n)/4^(k-1))){
        m_bar = m_bar - 1
        h = h/4
        #print(m_bar)
        return(list(est = est[m_bar,], m_bar = m_bar, h_m = h))
      }
    }
  }
  h = min(h,0.5)
  m_bar = M
  #print(m_bar)
  return(list(est = est[m_bar,], m_bar = m_bar, h_m = h))
}

FOMT <- function(Y, alpha = 0.05, h = NA, beta = NA, L = 1, sd = 0.1,R = 20, 
                 kernel = "epanechnikov",scale = 1, scale_A = 0.5){
  n = length(Y)
  if(!is.na(beta)){
    if((beta!=1) && (beta!=2)){
      print("Error! beta must be 1 or 2!")
      return(0)
    }
    if(is.na(h)){
      h = 0.5*(log(n)/n)^(1/(2*beta+1))
      #print(paste("Bandwidth h =",h))
    }
  }else{
    if(!is.na(h)){
      beta = 1
    }else{
      #tic()
      #print("Adaptive FOMT is computing...")
      C_s = ceiling(-2*log(alpha/2)*n/log(n))
      #P = matrix(NA,C_s*ceiling(R*log(n))*1*ceiling(log2(n/2)),2)
      P = matrix(NA,C_s*ceiling(R*log(n))*4*ceiling(log2(n/2)),2)
      number = 1
      A = rep(FALSE,n)
      #print("Generating indices set...")
      for (i in 1:C_s) {
        I = sample(1:n,1)
        if(I >= 2){
          for (m in 1:ceiling(R*log(n))){
            for (k_1 in 0:floor(log2(I-1))) {
              J = sample(1:2^k_1,1)
              A[I-J] = TRUE
              A[I] = TRUE
              P[number,1] = I-J
              P[number,2] = I
              number = number + 1
            }
          }
        }
        if(I <= n-1){
          for (m in 1:ceiling(R*log(n))) {
            for (k_2 in 0:floor(log2(n-I))) {
              J = sample(1:2^k_2,1)
              A[I] = TRUE
              A[I+J] = TRUE
              P[number,1] = I
              P[number,2] = I+J
              number = number + 1
            }
          }
        }
      }
      #print("Tuning bandwidth...")
      number = number - 1
      A = which(A == TRUE)
      #toc()
      #tic()
      res = CALM(Y,A,kernel,kappa = 1.5,sd = sd, scale_A = scale_A)
      Y_hat = res$est
      h = res$h_m
      #print(paste("m_bar = ", res$m_bar))
      #print(paste("Bandwidth by Lepskii h = ", res$h_m))
      temp = KL_max(kernel)
      K_max = temp[1]
      L_K = temp[2]
      lambda_0 = 0.5
      L_1 = sqrt(exp(1))*(K_max+L_K) 
      L_2 = (2*K_max+L_K)
      C_1 = (L_1/(8*L_2)+sqrt(exp(1))*K_max/lambda_0)^2
      C_star = 1
      W = max(C_1,C_star)
      N_max = number
      C_max = sd*sqrt(-2*log(alpha/N_max))*sqrt(C_star)/sqrt(n*h)
      #print("Searching potential violation...")
      #toc()
      #tic()
      for (i in 1:number) {
        I = P[i,1]
        J = P[i,2]
        I_Idx = which(A == I)
        J_Idx = which(A == J)
        est_1 = res$est[I_Idx]
        est_2 = res$est[J_Idx]
        dist = abs(J-I)/n*8*L_2/lambda_0
        if((dist <= h)&&(J/n<= (1-h))&&(I/n>=h)){
          critical_value = scale * min(C_max,C_max/h*dist/sqrt(C_star*C_1))
        }else{
          critical_value = scale*C_max
        }
        if((est_1-est_2) >= critical_value){
          # print("Violation between")
          # print(c(I,J))
          # print("found!")
          return(1)
        }
      }
      return(0)
      #tic()
    }
  }
  #print("beta")
  x = seq(1/n,1,1/n)

  C_s = ceiling(-2*log(alpha/2)/h)
  N_max = ceiling(40/log(2)*C_s*log(n)*log(n))
  
  if(beta == 1){
    lambda_0 = 0.5
  }
  if(beta == 2){
    lambda_0 = 0.03
  }
  temp = KL_max(kernel)
  K_max = temp[1]
  L_K = temp[2]
  L_1 = sqrt(exp(1))*(K_max+L_K) 
  L_2 = (2*K_max+L_K)*beta
  C_1 = (L_1/(8*L_2)+sqrt(exp(1))*K_max/lambda_0)^2
  C_star = 1
  W = max(C_1,C_star)

  C_max = sd*sqrt(-2*log(alpha/N_max))*sqrt(C_star)/sqrt(n*h)
  D_max = 16*K_max/lambda_0*L*h^beta
  middel_weight = weight_vector(0.5, n, h, beta, kernel)
  
  for (i in 1:C_s) {
    I = sample(1:n,1)
    if(I >= 2){
      for (m in 1:ceiling_FOMT(R*log(n))){
        for (k_1 in 0:floor_FOMT(log2(I-1))) {
          J = sample(1:2^k_1,1)
          est_1 = LP(I-J,Y,h,beta,kernel, middel_weight)
          est_2 = LP(I,Y,h,beta,kernel,middel_weight)
          dist = J/n*8*L_2/lambda_0
          #dist = J/n
          if(dist <= h){
            if((I/n<= 1-h)&&((I-J)/n>=h)){
              critical_value = scale * min(C_max,C_max/h*dist/sqrt(C_star*C_1))
              #critical_value = C_max/h*dist
            }else{
              critical_value = scale * min(C_max,C_max/h*dist/sqrt(C_star*C_1)) + D_max
              #critical_value = C_max/h*dist + D_max
            }
          }else{
            if((I/n<= 1-h)&&((I-J)/n>=h)){
              critical_value = scale * C_max
              #critical_value = C_max
            }else{
              critical_value = scale * C_max + D_max
              #critical_value = C_max + D_max
            }
          }
          if((est_1-est_2) >= critical_value){
            #print("Violation between")
            #print(c(I-J,I))
            #print("found!")
            return(1)
          }
        }
      }
    }
    if(I <= n-1){
      for (m in 1:ceiling_FOMT(R*log(n))) {
        for (k_2 in 0:floor_FOMT(log2(n-I))) {
          J = sample(1:2^k_2,1)
          est_1 = LP(I,Y,h,beta,kernel, middel_weight)
          est_2 = LP(I+J,Y,h,beta,kernel,middel_weight)
          dist = J/n*8*L_2/lambda_0
          #dist = J/n
          if(dist <= h){
            if(((I+J)/n<= 1-h)&&(I/n>=h)){
              critical_value = scale * min(C_max,C_max/h*dist/sqrt(C_star*C_1))
              #critical_value = C_max/h*dist
            }else{
              critical_value = scale * min(C_max,C_max/h*dist/sqrt(C_star*C_1)) + D_max
              #critical_value = C_max/h*dist + D_max
            }
          }else{
            if(((I+J)/n<= 1-h)&&(I/n>=h)){
              critical_value = scale * C_max
              #critical_value = C_max
            }else{
              critical_value = scale * C_max + D_max
              #critical_value = C_max + D_max
            }
          }
          if((est_1-est_2) >= critical_value){
            #print("Violation between")
            #print(c(I,I+J))
            #print("found!")
            return(1)
          }
        }
      }
    }
  }
  return(0)
}





