library(glmnet)

X=rtPCRgene_no6h
Y=mito_MP


test<-cv.glmnet(X,Y)

print(test)
plot(test$cvm)
summary(test)

cv1<-cv.glmnet(X,Y,alpha=1)
cv.5<-cv.glmnet(X,Y,alpha=0.5)
cv0<-cv.glmnet(X,Y,alpha=0.3)

par(mfrow=c(2,2))
plot(cv1,xlim=c(-6,0));plot(cv.5,xlim=c(-6,0));plot(cv0,xlim=c(-6,0))
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name,xlim = c(-3,0),ylim = c(0.01,0.03))
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey",ylim=c(0.01,0.03))
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue",ylim=c(0.01,0.03))
legend("topright",legend=c("alpha= 1","alpha= .5","alpha= .3"),pch=19,col=c("red","grey","blue"))
