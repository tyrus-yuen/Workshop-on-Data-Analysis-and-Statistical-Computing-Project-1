Gibbs_MH_xMatrix=read.csv(file="Gibbs_MH_xMatrix.csv") 
Gibbs_MH_yMatrix=read.csv(file="Gibbs_MH_yMatrix.csv")

dim(Gibbs_MH_yMatrix) 

x=Gibbs_MH_xMatrix[,2:101]
y=Gibbs_MH_yMatrix[,2:101]


chain=cbind(x[,1],y[,1])
colnames(chain)=c("x","y")
par(mfrow=c(2,1))
plotAuto(chain)


x_thin5=x[seq(from=1,to=99000,by=5),]
y_thin5=y[seq(from=1,to=99000,by=5),]
dim(x_thin5)


corrs_GMH_thin5=rep(NA,100)
for(i in 1:100){
  corrs_GMH_thin5[i]=cor(x_thin5[,i],y_thin5[,i])
  }


mean(corrs_GMH_thin5)
sd(corrs_GMH_thin5)
LCI_GMH_thin5=mean(corrs_GMH_thin5)-qnorm(0.975,mean=0,sd=1)*sd(corrs_GMH_thin5)
UCI_GMH_thin5=mean(corrs_GMH_thin5)+qnorm(0.975,mean=0,sd=1)*sd(corrs_GMH_thin5)

par(mfrow=c(1,1))
plot(corrs_GMH_thin5,xlab="MEAN=-0.06070848 SD=0.005857476")
abline(h=mean(corrs_GMH_thin5),col="red")
abline(h=LCI_GMH_thin5,col="blue")
abline(h=UCI_GMH_thin5,col="blue")



chain_thin5=cbind(x[,1],y[,1])
colnames(chain_thin5)=c("x","y")

plotAuto(chain)
plotAuto(chain, thin=5)
plotAuto(chain_thin5)


write.csv(x_thin5,file="Gibbs_MH_thin5_xMatrix.csv")
write.csv(y_thin5,file="Gibbs_MH_thin5_yMatrix.csv")