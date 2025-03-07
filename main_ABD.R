source("Phi_ABD.R")
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


library(tictoc)
alpha = 0.05
M = 100
sd = 0.3
#n = c(100,200,500,1000,2000,5000,10000)
n = c(400,800,1200,1600,2000,2400,2800,3200)
cv = c(0.1109668,0.07742796,6.675786e-02,5.876668e-02,5.388069e-02,
       4.80369e-02,4.37757e-02,4.203237e-02) 
# critical values for n = 400, 800, 1200, 1600, 2000, 2400, 2800 and 3200.
R = 100
for (m in 1:8) {
  x = seq(1/n[m],1,1/n[m])
  if(!is.na(cv[m])){
    c = cv[m]
    print(c(n[m],c))
  }else{
    tic()
    c = Simulate_c(n[m],alpha,R)
    toc()
    print(c(n[m],c))
  }
  for (i in c(0,3,5,4,15)) {
    # Test functions 0, 3, 5, 4, 15 correspond to f0,f1,f2,f3 and f4 in our paper,
    # respectively
    signal = test_functions(x,i,sd)
    res = 0
    tic()
    tm = microbenchmark(rejection = {res = res + Phi_Durot(Y = signal + rnorm(n[m],0,sd),
                                                           alpha = alpha, sd = sd,c = c)},
                        times = M, unit = "s")
    toc()
    saveRDS(tm,paste("Time of F",i," with size ",n[m]))
    print(paste("Rejection of F",i," with size ",n[m], " = ", res))
  }
}

