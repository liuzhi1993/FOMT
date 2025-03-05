source("FOMT.R")
test_functions = function(x,i=0,sd=0.1,b=0.15){
  n = length(x)
  if(i==0){
    res = rep(0,n)
  }
  if(i==1){
    res = 1 + x - 0.15*exp(-50*(x-0.5)^2)
  }
  if(i==2){
    res = 1 + x - 0.25*exp(-50*(x-0.5)^2)
  }
  if(i==3){
    res = 1 + x - 0.45*exp(-50*(x-0.5)^2)
  }
  if(i==4){
    #res = -1.5*sd*x
    res = -1.5*0.2*x
  }
  if(i==5){
    res = -0.2*exp(-50*(x-0.5)^2)
  }
  if(i==6){
    res = 0.1*cos(6*pi*x)
  }
  if(i==7){
    res = 0.2*x+test_functions(x,i=5)
  }
  if(i==8){
    res = 0.2*x+test_functions(x,i=6)
  }
  if(i==9){
    res = x*(1-x)
  }
  if(i==10){
    res = exp(-100*(x-0.25)^2)-10*(x-0.5)^3*(x<=0.5)-0.1*(x-0.5)*(x>0.5)
  }
  if(i==11){
    res = 0.3*(x-0.5)-exp(-250*(x-0.25)^2)+15*(x-0.5)^3*(x<=0.5)
  }
  if(i==12){
    res = x + 0.416*exp(-50*x^2)
  }
  if(i==13){
    res = x-4*pnorm(1*x)
  }
  if(i==14){
    res = x-1.2*pnorm(5*x)
  }
  if(i==15){
    res = x -1.5*pnorm(4*x)
  }
  return(res)
}

# FOMT
library(microbenchmark)
library(tictoc)
alpha = 0.05
M = 1000
sd = 0.3
n = c(400,800,1200,1600,2000,2400,2800,3200)
R = 0.1
scale = 0.58
for (m in 1:length(n)) {
  x = seq(1/n[m],1,1/n[m])
  for (i in c(0,3,5,4,9)) { 
    # Test functions 0, 3, 5, 4, 9 correspond to f0,f1,f2,f3 and f4 in our paper,
    # respectively
    signal = test_functions(x,i,sd)
    res = 0
    tic()
    tm = microbenchmark(rejection = {res = res + FOMT(Y = signal+rnorm(n[m],0,sd), 
                                                      alpha = alpha, beta = 2, sd = sd, R = R, scale = scale)},
                        times = M, unit = "s")
    toc()
    #saveRDS(tm,paste("New Time of F",i," with size ",n[m],"scale = ", scale, "R =",R))
    print(paste("Rejection of F",i," with size ",n[m], " = ", res))
  }
}

R = 0.1
scale = 1
scale_A = 0.00175
for (m in 8:length(n)) {
  x = seq(1/n[m],1,1/n[m])
  for (i in c(0,3,5,4,9)) {
    signal = test_functions(x,i,sd)
    print(max(abs(diff(signal)*n[m])))
    res = 0
    tic()
    tm = microbenchmark(rejection = {res = res + FOMT(Y = signal+rnorm(n[m],0,sd), 
                                                      alpha = alpha, beta = NA, 
                                                      sd = sd, R = R, scale = scale,
                                                      scale_A = scale_A, power = 1.6,
                                                      adaptivity = T)},
                        times = M, unit = "s")
    toc()
    #saveRDS(tm,paste("New LP Time of F",i," with size ",n[m],"scale = ", scale,"scale_A = ", scale_A, "R =", R))
    print(paste("LP Rejection of F",i," with size ",n[m], " = ", res))
  }
}

