#Prepares dvh data for fitting
prepare.DVH<-function(DVH,n=NULL)
{
  len<-length(DVH)
  d<-matrix(rep(0,1000*len),ncol = len)
  v<-matrix(rep(0,1000*len),nrow=len)
  binSize<-rep(NA,len)
  dvh<-lapply(DVH, function(x)
    {
    x[,2]<-1-(x[,2]/max(x[,2]))
    x[,2]<-ifelse(x[,2]==0,1e-10,x[,2])
    x[,2]<-ifelse(x[,2]==1,1-1e-10,x[,2])
    x[,1]<-stats::qnorm(x[,2])
    x[,2]<-stats::dnorm(x[,1])
    x
  })

  if(is.null(n))
  {
    len2<-30
    n<-seq(0.1,1,length.out = len2)
  }
  else
  {
    if(sum(n<0 | n>1)>0)stop("n must lie within 0 and 1")
    len2<-length(n)
  }

#Normal transformation so we get uncorrelated odd moments

  EUDn<-matrix(rep(NA,len2*len),ncol=len2)
  for(k in 1:len2)
  {

    EUDn[,k]<- unlist(lapply(dvh, function(x,k)(sum(((x[,1]+10)^(1/n[k]))*x[,2]))^n[k],k=k))
  }
maxEUD<-max(EUDn)
  list(EUD=EUDn/maxEUD
       ,maxEUD=maxEUD, n=n)
}


