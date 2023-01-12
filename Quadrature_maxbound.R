#Quadrature
f3 <- function(x) {return((abs(cos(x))/x)*exp(-(log(x)-3)^2)) }#function of question 3
f3_1 <- function(x) {}
body(f3_1) <- D(body(f3), 'x')

f3_2 <- function(x) {}
body(f3_2) <- D(body(f3_1), 'x')

f3_3 <- function(x) {}
body(f3_3) <- D(body(f3_2), 'x')

f3_4 <- function(x) {}
body(f3_4) <- D(body(f3_3), 'x') #fourth derivative
curve(f3, 0, 200 , ylab = "y=f(x)") #sketch the curve
curve(f3, 0.0005, 5000 , ylab = "y=f(x)") #sketch the curve

#1.Composite trapezoidal rule
CTrape<-function(funct,a,b,n){ #Define the parameters #n is the no. of subinterval
  h <- (b - a) / n #Calculate the increment
  c <- 1:(n - 1) #Create a vector for the order of the increment 
  xj <- a + c * h #Transform the order vector to a vector of f(x_j)
  Intergrand <- (1/2*h) * (2 * sum(funct(xj))) 
  #Composite Trapezoidal Rule= h/2+[f(a)+Summation_c=1_n-1(f(x_j)+f(b))]
  
  value <- function (funct, x, adj=h/1000000) { #Set up an adjustment valueCompromise on a small error to avoid the infinity value on 0
    if (is.nan(funct(x))) { #if condition when f(x)=NaN
      if (x==0) return (funct(adj+x)) #add adjustment value to x_1 when x=0
      else return (func(adj*(-1)+x)) #prevent any NaN cases of x_j!=0
    }
    else {return(funct(x))}
  } 
  
  Intergrand=Intergrand+(h/2*(value(funct,a)+value(funct,b))) #Calculate the final intergrand

  (Result= Intergrand)
}
set.seed(77960)
CTrape(f3, 1/1e7, 1e7, 5000)

#2. Composite Simpson's rule
CSimp<-function(funct,a,b,n) { #Define the parameters #n is the no. of subinterval
  h <- (b - a) / n #Calculate the increment
  c<-0:(n-1)
  d <- 1:(n - 1) #Create a vector for the order of the increment 
  xj1 <- a+(c+1/2)*h  #Transform the order vector to a vector of f(x_j)
  xj2<- a+(d*h)
  Intergrand <- (1/6*h) * (4 * sum(funct(xj1))+2*sum(funct(xj2))) 
  #Calculate parts of the formula of Simpson's rule
  
  value <- function (funct, x, adj=h/1000000) { #Set up an adjustment valueCompromise on a small error to avoid the infinity value on 0
    if (is.nan(funct(x))) {
      if (x==0) return (funct(adj+x))
      else return (func(adj*(-1)+x))
    }
    else {return(funct(x))}
  } 
  
  Intergrand=Intergrand+(h/6*(value(funct,a)+value(funct,b)))
  #Composite Simpson's Rule= h/6*[f(a)+4*Summation_i=0_n-1(f(a+(i+1/2)*h))+2*Summation_i=1_n-1(f(a+i*h)+f(b)]
  return(Intergrand)
}
set.seed(77960)
CSimp(f3, 1/1e7, 1e7, 5000)


#Calculate the running time
runtime=function(fn){ 
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(-1*(start_time-end_time))
}
runtime(CTrape(f3, 1/1e7, 1e7, 5000))
runtime(CSimp(f3, 1/1e7, 1e7, 5000))

x2<-max(f3_2())
x4<-(f3_4())

((1e7-1/1e7)^3)/(12*(2000)^2)*x2 #Error Estimation of CTrape
((1e7-1/1e7)^5)/(180*(2000)^4)*x4 #Error Estimation of CSimp