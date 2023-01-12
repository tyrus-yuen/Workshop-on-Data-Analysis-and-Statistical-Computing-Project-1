MH_xMatrix=read.csv(file="MH_xMatrix.csv") 
MH_yMatrix=read.csv(file="MH_yMatrix.csv")

dim(MH_yMatrix) 

x=MH_xMatrix[,2:101]
y=MH_yMatrix[,2:101]


chain=cbind(x[,1],y[,1])
colnames(chain)=c("x","y")
par(mfrow=c(2,1))
plotAuto(chain)


x_thin5=x[seq(from=1,to=99000,by=5),]
y_thin5=y[seq(from=1,to=99000,by=5),]
dim(x_thin5)


corrs_MH_thin5=rep(NA,100)
for(i in 1:100){
  corrs_MH_thin5[i]=cor(x_thin5[,i],y_thin5[,i])
  }


mean(corrs_MH_thin5)
sd(corrs_MH_thin5)
LCI_MH_thin5=mean(corrs_MH_thin5)-qnorm(0.975,mean=0,sd=1)*sd(corrs_MH_thin5)
UCI_MH_thin5=mean(corrs_MH_thin5)+qnorm(0.975,mean=0,sd=1)*sd(corrs_MH_thin5)

par(mfrow=c(1,1))
plot(corrs_MH_thin5,xlab="MEAN=-0.06070848 SD=0.005857476")
abline(h=mean(corrs_MH_thin5),col="red")
abline(h=LCI_MH_thin5,col="blue")
abline(h=UCI_MH_thin5,col="blue")



chain_thin5=cbind(x[,1],y[,1])
colnames(chain_thin5)=c("x","y")

plotAuto(chain)
plotAuto(chain, thin=5)
plotAuto(chain_thin5)

write.csv(corrs_MH_thin5,file="correlation_MH_thin5.csv")
write.csv(x_thin5,file="MH_thin5_xMatrix.csv")
write.csv(y_thin5,file="MH_thin5_yMatrix.csv")