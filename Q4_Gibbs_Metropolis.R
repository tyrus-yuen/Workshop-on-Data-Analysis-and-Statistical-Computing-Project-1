install.packages("plotMCMC")
install.packages("MASS")
install.packages("plot3D")
library(plotMCMC)
library(MASS)
library(plot3D)


#install.packages("xlsx")
#install.packages("rJava")
#install.packages("xlsxjars")
#library(rJava)
#library(xlsxjars)
#library(xlsx)

runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}

# Gibbs + Metropolis

f=function(x,y){
  #if(x<0|x>1|y<0|y>1) return("x,y out of range")
  result=(sin(pi*x)*sin(pi*x^2)^20+sin(pi*y)*sin(2*pi*y^2)^20)
  return(result)
}

fx=function(x){
  #if(x<0|x>1) return("x,y out of range")
  result=sin(pi*x)*sin(pi*x^2)^20
  return(result)
}

fxx=function(x){
  #if(x<0|x>1) return("x,y out of range")
  result=pi*cos(pi*x)*sin(pi*x^2)^20+20*sin(pi*x^2)^19*cos(pi*x^2)*2*pi*x
  return(result)
}

fxy=function(x,y0){
  #if(x<0|x>1) return("x,y out of range")
  result=(sin(pi*x)*sin(pi*x^2)^20+fy(y0))/(Itgfx+fy(y0))
  return(result)
}


fy=function(y){
  #if(y<0|y>1) return("x,y out of range")
  result=sin(pi*y)*sin(2*pi*y^2)^20
  return(result)
}

fyy=function(y){
  #if(y<0|y>1) return("x,y out of range")
  result=pi*cos(pi*y)*sin(2*pi*y^2)^20+20*sin(2*pi*y^2)^19*cos(2*pi*y^2)*4*pi*y
  return(result)
}

fyx=function(y,x0){
  #if(y<0|y>1) return("x,y out of range")
  result=(sin(pi*y)*sin(2*pi*y^2)^20+fx(x0))/(Itgfy+fx(x0))
  return(result)
}

#Inverse CDF Method, to find the integral of fx and fy
n=1e7
ux=runif(n,0,1)
uy=runif(n,0,1)
Itgfx=mean(fx(ux))
Itgfy=mean(fy(uy))

#find maxfx and maxfy
#Biseciton Method
bisection=function(fx,a,b,error,step){
  i=0
  if(fx(a)*fx(b)>0) return("Invalid Interval [a,b]")
  
  while(b-a>error){
    c=(a+b)/2
    ifelse(fx(a)*fx(c)<0,b<-c,a<-c);
    i=i+1
    if(i>=step) return("Steps has been used up")
  }
  root=(a+b)/2
  return(root)
}

max_x=fx(bisection(fxx,0.2,0.8,1e-10,100))
max_y=fy(bisection(fyy,0,0.8,1e-10,100))


#Gibbs + Metropolis
#One dimentional Metropolis
Metropolis_1D=function(f,p1,x0,n){
  x=rep(NA,n)
  i=1
  x_c=x0
  while(i<=n){
    x_c=p1(x0)
    u=runif(1,0,1)
    a=f(x_c)/f(x0)
    if(u<=a){
      x[i]=x_c
      }
    else{
      x[i]=x0
    }
    x0=x[i]
    i=i+1
  }
  chain=x
  result=list("chain"=chain)
  return(result)
}


p1=function(x){
  result=runif(1,0,1)
  return(result)
}

Gibbs_MH=function(fx,fy,x0,y0,n){
  x=rep(NA,n)
  y=rep(NA,n)
  for (i in 1:n){
    fyx=function(y){
      #if(y<0|y>1) return("x,y out of range")
      result=(sin(pi*y)*sin(2*pi*y^2)^20+fx(x0))/(Itgfy+fx(x0))
      return(result)
    }
    #mx=(fx(x0)+max_y)/(fx(x0)+Itgfy)
    chainy=Metropolis_1D(fyx,p1,y0,1)
    y[i]=chainy$chain
    y0=y[i]
    fxy=function(x){
      #if(x<0|x>1) return("x,y out of range")
      result=(sin(pi*x)*sin(pi*x^2)^20+fy(y0))/(Itgfx+fy(y0))
      return(result)
    }
    #mx=(max_x+fy(y0))/(Itgfx+fy(y0))
    chainx=Metropolis_1D(fxy,p1,x0,1)
    x[i]=chainx$chain
    x0=x[i]
  }
  chain=cbind(x,y)
  result=list("chain"=chain)
  return(result)
}

chainxy=Gibbs_MH(fx,fy,0.5,0.5,1e5)
x=chainxy$chain[(1+1e3):1e5,1]
y=chainxy$chain[(1+1e3):1e5,2]
cor(x,y)

z=rep(NA,length(x))
for(i in 1:length(x)){
  z[i]=f(x[i],y[i])
}
scatter3D(x,y,z,colvar=z,pch=19,cex=0.5)
plotAuto(chainxy$chain)
plotAuto(chainxy$chain,thin=5)


#repeated trails for m=100 times
m=100
n=1e5
burn=1e3
corrs_GMH=rep(NA,m)
xmatrix_GMH=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)
ymatrix_GMH=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)
starttime=Sys.time()
for (i in 1:m){
  chainxy=Gibbs_MH(fx,fy,0.5,0.5,n)
  x=chainxy$chain[(1+burn):n,1]
  y=chainxy$chain[(1+burn):n,2]
  xmatrix_GMH[,i]=x
  ymatrix_GMH[,i]=y
  corrs_GMH[i]=cor(x,y)
}
endtime=Sys.time()
runtime_GMH_100=endtime-starttime
runtime_GMH_100
runtime_GMH=runtime(Gibbs_MH(fx,fy,0.5,0.5,n))
runtime_GMH

mean(corrs_GMH)
sd(corrs_GMH)
plot(corrs_GMH)
range(corrs_GMH)
max(corrs_GMH)-min(corrs_GMH)

LCI_GMH=mean(corrs_GMH)-qnorm(0.975,mean=0,sd=1)*sd(corrs_GMH)
UCI_GMH=mean(corrs_GMH)+qnorm(0.975,mean=0,sd=1)*sd(corrs_GMH)
LCI_GMH
UCI_GMH
LenCI_GMH=UCI_GMH-LCI_GMH
LenCI_GMH

write.csv(xmatrix_GMH,file="Gibbs_MH_xMatrix.csv")
write.csv(ymatrix_GMH,file="Gibbs_MH_yMatrix.csv")
write.csv(corrs_GMH,file="Gibbs_MH_correlation.csv")
write.csv(c(runtime_GMH,runtime_GMH_100),file="Gibbs_MH_runtime.csv")