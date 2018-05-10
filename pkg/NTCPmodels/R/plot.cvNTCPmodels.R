plot.cvNTCPmodels<-function(object,...)
{
  xx<-as.numeric(object[[1]][,1])
  yy<-as.numeric(object[[1]][,2])
  if(object$type.measure=="auc")ylab="AUC"
  graphics::plot(object[[1]],
       xlab = "n",cex.lab=1.3,type="b",
       ylab = ylab,ylim=range(c(yy+stats::qnorm(0.972)*object$sd,
                                yy-stats::qnorm(0.972)*object$sd)),
       pch=16,col="red")
  graphics::arrows(xx,yy-stats::qnorm(0.972)*object$sd,xx,
                      yy+stats::qnorm(0.972)*object$sd,
                   length = 0.08,angle = 90,code = 3,col="gray")
  graphics::abline(v=xx[yy==max(yy)],lty=2)
}


