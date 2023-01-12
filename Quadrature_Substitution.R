#Quadrature by substitution

fy=function(y){
  result=(1/(y*(1-y)))*abs(cos((1-y)/y))*exp(-(log((1-y)/y)-3)^2)
  return(result)
}
curve(fy,0,1, ylab = "y=f(x)") #sketch the curve in 0 to 1 
curve(fy,0.0005,1/2000*10, ylab = "y=f(x)") #sketch the curve in 0 to 0.1 

#1.Composite trapezoidal rule
CTrape<-function(funct,a,b,n){ #Define the parameters #n is the no. of subinterval
  h <- (b - a) / n #Calculate the increment
  c <- 1:(n - 1) #Create a vector for the order of the increment 
  xj <- a + c * h #Transform the order vector to a vector of f(x_j)
  Integral <- (1/2*h) * (2 * sum(funct(xj))) 
  #Composite Trapezoidal Rule= h/2+[f(a)+Summation_c=1_n-1(f(x_j)+f(b))]
  
  value<-function (funct, x, adj=h/1e6) { #Set up an adjustment valueCompromise on a small error to avoid the infinity value on 0
    if (is.nan(funct(x))) { #if condition when f(x)=NaN
      if (x==0) return (funct(adj+x)) #add adjustment value to x_1 when x=0
      else return (funct(adj*(-1)+x)) #prevent any NaN cases of x_j!=0
    }
    else {return(funct(x))}
  } 
  
  Integral=Integral+(h/2*(value(funct,a)+value(funct,b))) #Calculate the final Integral
  (Result= Integral)
}
set.seed(77960)
CTrape(fy,0.0005,1,2000)

#2. Composite Simpson's rule
CSimp<-function(funct,a,b,n) { #Define the parameters #n is the no. of subinterval
  h <- (b - a) / n #Calculate the increment
  c<-0:(n-1)
  d <- 1:(n - 1) #Create a vector for the order of the increment 
  xj1 <- a+(c+1/2)*h  #Transform the order vector to a vector of f(x_j)
  xj2<- a+(d*h)
  Integral <- (1/6*h) * (4 * sum(funct(xj1))+2*sum(funct(xj2))) 
  #Calculate parts of the formula of Simpson's rule
  
  value <- function (funct, x, adj=h/1000000) { #Set up an adjustment valueCompromise on a small error to avoid the infinity value on 0
    if (is.nan(funct(x))) {
      if (x==0) return (funct(adj+x))
      else return (funct(adj*(-1)+x))
    }
    else {return(funct(x))}
  } 
  
  Integral=Integral+(h/6*(value(funct,a)+value(funct,b)))
  #Composite Simpson's Rule= h/6*[f(a)+4*Summation_i=0_n-1(f(a+(i+1/2)*h))+2*Summation_i=1_n-1(f(a+i*h)+f(b)]
  (Result= Integral) #return the result
}
set.seed(77960)
CSimp(fy,0.0005,1,2000)

#Calculate the running time
runtime=function(fn){ 
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}


runtime(CTrape(fy,5e-4,1,2000))

runtime(CSimp(fy,5e-4,1,2000))


