plot.cvNTCPmodels<-function(object)
{
  plot(as.numeric(unlist(object[1])),as.numeric(unlist(object[2])),xlab = "n",ylab = names(object)[2],pch=16,col="red")
}

