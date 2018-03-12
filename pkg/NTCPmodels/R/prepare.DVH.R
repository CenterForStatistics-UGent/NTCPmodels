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
    x[,1]<-qnorm(x[,2])
    x[,2]<-dnorm(x[,1])
    x
  })
 #  par(mfrow=c(1,2))
 # plot(dvh[[200]][,1],dvh[[200]][,2],type="l")
 # lines(dvh[[6]][,1],dvh[[6]][,2],type="l")
 #
 # plot(DVH[[200]][,1], c(0,diff(1-(DVH[[200]][,2]/max(DVH[[200]][,2])))),type="l")
 # lines(DVH[[6]][,1],c(0,diff(1-(DVH[[6]][,2]/max(DVH[[6]][,2])))),type="l")
 #  for(i in 1:len)
 #  {
 #    dvh<-DVH[[i]]
 #    lenn<-length(dvh[,2])
 #    v[i,][1:(lenn)]<-c(0,-diff(dvh[,2]/max(dvh[,2])))
 #    grp<-diff(dvh[,1])
 #    if(unique(grp)>1)stop(paste("The dose bin size of the observation",i,"has a problem"))
 #  }
  # binSize<-as.factor(binSize)
  # levs<-as.numeric(levels(binSize))
  # ll1<-length(levs)
  # doseBins<-matrix(rep(0,ll1*1000),ncol=ll1)
  # for(ii in 1:ll1)
  # {
  #   seqq<-seq(levs[ii],100,by=levs[ii])
  #   doseBins[,ii][1:length(seqq)]<-seqq
  # }
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

