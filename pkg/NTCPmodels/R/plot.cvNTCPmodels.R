plot.cvNTCPmodels<-function(object)
{
  xx<-as.numeric(object[[1]][,1])
  yy<-as.numeric(object[[1]][,2])
  if(object$type.measure=="auc")ylab="AUC"
  plot(object[[1]],
       xlab = "n",cex.lab=1.3,type="b",
       ylab = ylab,ylim=range(c(yy+qnorm(0.972)*object$sd,yy-qnorm(0.972)*object$sd)),
       pch=16,col="red")
  arrows(xx,yy-qnorm(0.972)*object$sd,xx,yy+qnorm(0.972)*object$sd,length = 0.08,angle = 90,code = 3,col="gray")
  abline(v=xx[yy==max(yy)],lty=2)
}


