source("Phi_DS.R")
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
    res = -1.5*sd*x
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

# DS
library(microbenchmark)
library(tictoc)
alpha = 0.05
M = 100
sd = 0.3
n = c(400,800,1200,1600,2000,2400,2800,3200)
R = 100
for (m in 1:8) {
  x = seq(1/n[m],1,1/n[m])
  print(n[m])
  tic()
  critical_value = Simulate_c(n = n[m], R = 100, type = 1)
  print(c(n[m],critical_value))
  toc()
  for (i in c(3,7,4,15,0)) {
    signal = -test_functions(x,i,sd)
    res = 0 
    tic()
    tm = microbenchmark(rejection = {res = res + Phi_DS(Y = signal+rnorm(n[m],0,sd), 
                                     sd = sd, type = 1, critical_value = critical_value)},
                          times = M, unit = "s")
    toc()
    #saveRDS(tm,paste("DS 1 Time of F",i," with size ",n[m]))
    print(paste("DS 1 Rejection of F",i," with size ",n[m], " = ", res))
  }
}

for (m in 1:6) {
  x = seq(1/n[m],1,1/n[m])
  # tm_simulation = microbenchmark(Simulate_c(n = n[m],R = 100, type = 2),
  #                                times = 10, unit = "s")
  # saveRDS(tm_simulation,paste("2 Time of simulation with size ",n[m]))
  print(n[m])
  tic()
  critical_value = Simulate_c(n = n[m], R = 100, type = 2)
  toc()
  for (i in 0:15) {
    signal = -test_functions(x,i,sd)
    res = 0 
    tm = microbenchmark(rejection = {res = res + Phi_DS(Y = signal+rnorm(n[m],0,sd), 
                                                        sd = sd, type = 2, critical_value = critical_value)},
                        times = M, unit = "s")
    saveRDS(tm,paste("DS 2 Time of F",i," with size ",n[m]))
    print(paste("DS 2 Rejection of F",i," with size ",n[m], " = ", res))
  }
}



