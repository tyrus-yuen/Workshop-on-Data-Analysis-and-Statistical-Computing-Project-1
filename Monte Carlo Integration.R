#Q3 Monte Carlo Integration & Quadrature
#I=Integrate_{0}^{Infinity} (abs|cos(x)|/x*exp^(-(log(x)-3)^2)) dx
#f(x)=exp^(-x) g(x)/f(x)= Integrate_{0}^{Infinity} (abs|cos(x)|/x*exp^(-(log(x)-3)^2+x)) dx
f3<- function(X){abs(cos(X))/X*exp(-1*((log(X)-3)^2))} #Original Formula
f3a <- function(X){abs(cos(X)/X*exp(-1*(log(X)-3)^2)/exp(-X))} #Pick Random Variable X following Exp~(1)
f3b <- function(X){abs(cos(X))*sqrt(2*pi)*exp(-1/2*((log(X)-3)^2))} #Pick Random Variable X following Lognormal~(3,1)
curve(f3, 0, 200,xlab= "Q3 function", ylab = "y=f(x)") #plot the curve
#g(x)/f(x)=abs|cos(x)|/x*exp^(-(log(x)-3)^2+x)

#Monte Carlo Integration
#Exponential
mcie <- function(n,f){          #Set up a function of simulation
  U <- runif(n,0,1)    # Simulate n random numbers 
  X <- -log(U)/1 #Simulate n values from Exp(1)
  Int <- c(mean(f(X)), var(f(X))/n) #Compute g(x_i)/f(x_i)=x_1 for i = 1 to n and the error estimation
  Int} #Calculate the sample mean to estimate the population mean
set.seed(77960)
mcie(5000000,f3a)

#Lognormal 
mciln <- function(n,f){          #Set up a function of simulation
  X <- rlnorm(n, meanlog=3, sdlog=1)  #Simulate n values from Lognormal(0,1)
  Int <- c(mean(f(X)), var(f(X))/n) #Compute g(x_i)/f(x_i)=x_1 for i = 1 to n and the error estimation
  Int} #Calculate the sample mean to estimate the population mean
set.seed(77960)
mciln(5000000,f3b)


#Calculate the running time
runtime=function(fn){ 
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(start_time-end_time)
}
runtime(mcie(100000,f3))
runtime(mciln(100000,f3))




