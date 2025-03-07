source("Phi_C.R")
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

# C
library(microbenchmark)
library(tictoc)
alpha = 0.05
M = 100
sd = 0.3
n = c(400,800,1200,1600,2000,2400,2800,3200)
cv = c(3.07508,3.501234,3.279791,3.293203,3.222466,3.402945,3.490897,2.994733) 
# critical values for n = 400, 800, 1200, 1600, 2000, 2400, 2800 and 3200.
R = 100
kernel = "epanechnikov"
k = 0
for (m in 1:8) {
  x = runif(n[m])
  x = seq(1/n[m],1,1/n[m])
  #tic()
  #c = simulate_c(n = n[m],X = x, R = 20, kernel = kernel, k=k, alpha = alpha)
  #print(c(n[m],c))
  #toc()
  c = cv[m]
  for (i in c(0,3,5,4,15)) {
    # Test functions 0, 3, 5, 4, 15 correspond to f0,f1,f2,f3 and f4 in our paper,
    # respectively
    x = runif(n[m])
    signal = test_functions(x,i,sd)
    res = 0 
    tic()
    tm = microbenchmark(rejection = {res = res + Phi_C(X = x, Y = signal+rnorm(n[m],0,sd),
                                                       kernel="epanechnikov",sd = sd, k = k,c = c)},
                        times = M, unit = "s")
    toc()
    saveRDS(tm,paste("Time of F",i," with size ",n[m]))
    print(paste("Rejection of F",i," with size ",n[m], " = ", res))
  }
}
