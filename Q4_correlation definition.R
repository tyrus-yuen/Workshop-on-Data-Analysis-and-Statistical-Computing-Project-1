runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}


#Monte Carlo
MC <-function(f,n){
 x<-rep(NA,n);
 i=1;
 r <- runif(n);
 while(i<=n){
 x[i] <- f(r[i]);
 i=i+1;
 }
 result=mean(x);
 return(result);
}
  
  #let f(x,y)= f(x) + f(y)
  fx<-function(x){sin(x*pi)*sin(pi*x^2)^20}
  fy<-function(y){sin(y*pi)*sin(2*pi*y^2)^20}
  xfx<-function(x){x*sin(x*pi)*sin(pi*x^2)^20}
  yfy<-function(y){y*sin(y*pi)*sin(2*pi*y^2)^20}
  xxfx<-function(x){(x^2)*sin(x*pi)*sin(pi*x^2)^20}
  yyfy<-function(y){(y^2)*sin(y*pi)*sin(2*pi*y^2)^20}

  1e8

Corr_Q4=function(N){  
  start_time=Sys.time()
  mfx <- MC(fx,N)
  mxfx <- MC(xfx,N)
  mxxfx <- MC(xxfx,N)
  mfy <- MC(fy,N)
  myfy <- MC(yfy,N)
  myyfy <- MC(yyfy,N)
  
  c <- mfx + mfy
  
  #E(XY) = 1/2F(xf(x)) + 1/2F(yf(y))
  #E(X) = F(xf(x)) + 1/2F(f(y))
  #E(X^2) = F(x^2f(x)) + 1/3F(f(y))
  #E(Y) = F(yf(y)) + 1/2F(f(x))
  #E(Y^2) = F(y^2f(y)) + 1/3F(f(x))
  
  Ex<-(mxfx + mfy/2)/c
  Ey<-(myfy + mfx/2)/c
  Exx<-(mxxfx + mfy/3)/c
  Eyy<-(myyfy + mfx/3)/c
  Exy<-(myfy/2 + mxfx/2)/c
  
  corr<-(Exy-Ex*Ey)/sqrt((Exx-Ex^2)*(Eyy-Ey^2))
  end_time=Sys.time()
  runtime=end_time-start_time
  result=list("corr"=corr,"runtime"=runtime)
  return(result)
}


m=100
corr_def=rep(NA,m)
corr_runtime=rep(NA,m)
for(i in 1:m){
result=Corr_Q4(1e7)
corr_def[i]=result$corr
corr_runtime[i]=result$runtime
}
corrs_def=corr_def

mean(corrs_def)
sd(corrs_def)
LCI_def=mean(corrs_def)-qnorm(0.975,mean=0,sd=1)*sd(corrs_def)
UCI_def=mean(corrs_def)+qnorm(0.975,mean=0,sd=1)*sd(corrs_def)

plot(corr_def,xlab="MEAN=-0.06101332 SD=0.0007891465")
abline(h=mean(corrs_def),col="red")
abline(h=LCI_def,col="blue")
abline(h=UCI_def,col="blue")

write.csv(corrs_MH,file="MH_correlation_thin.csv")
write.csv(c(runtime_MH,runtime_MH_100),file="MH_runtime_thin.csv")



#Correlation=NULL
#plot(Correlation,xlab="",xlim=c(0,100),ylim=c(-0.0598,-0.0611))
#abline(h=corr_ture,col="green")
#abline(h=mean(corrs_MH),col="red")
#abline(h=-0.05990,col="yellow")
#abline(h=-0.06023,col="blue")
#abline(h=-0.06105,col="orange")
#abline(h=-0.06042,col="black")
#abline(h=-0.06071,col="purple")


# [1] -0.06118968
  

#Directly using package cubature
  
  library(cubature)
  f<-function(x){(sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20)}
  fx<-function(x){x[1]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
  fy<-function(x){x[2]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
  fxx<-function(x){(x[1]^2)*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
  fyy<-function(x){(x[2]^2)*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
  fxy<-function(x){x[2]*x[1]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
  
  c<- as.numeric(hcubature(f,c(0,0),c(1,1))['integral'])
  
  Ex<- as.numeric(hcubature(fx,c(0,0),c(1,1))['integral'])/c
  Ey<- as.numeric(hcubature(fy,c(0,0),c(1,1))['integral'])/c
  Exx<- as.numeric(hcubature(fxx,c(0,0),c(1,1))['integral'])/c
  Eyy<- as.numeric(hcubature(fyy,c(0,0),c(1,1))['integral'])/c
  Exy<- as.numeric(hcubature(fxy,c(0,0),c(1,1))['integral'])/c
  
  corr_ture<-(Exy-Ex*Ey)/sqrt((Exx-Ex^2)*(Eyy-Ey^2))
  corr_ture
#[1] -0.06095707
  
