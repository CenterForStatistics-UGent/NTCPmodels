#Prepares dvh data for fitting
prepare.DVH<-function(DVH,n=NULL)
{
  len<-length(DVH)
  d<-matrix(rep(0,1000*len),ncol = len)
  v<-matrix(rep(0,1000*len),nrow=len)
  binSize<-rep(NA,len)
  for(i in 1:len)
  {
    dvh<-DVH[[i]]
    lenn<-length(dvh[,2])
    v[i,][1:(lenn)]<-c(0,-diff(dvh[,2]/max(dvh[,2])))
    grp<-diff(dvh[,1])
    if(unique(grp)>1)stop(paste("The dose bin size of the observation",i,"has a problem"))
    binSize[i]<-unique(grp)
  }
  binSize<-as.factor(binSize)
  levs<-as.numeric(levels(binSize))
  ll1<-length(levs)
  doseBins<-matrix(rep(0,ll1*1000),ncol=ll1)
  for(ii in 1:ll1)
  {
    seqq<-seq(levs[ii],100,by=levs[ii])
    doseBins[,ii][1:length(seqq)]<-seqq
  }
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


  EUDn<-matrix(rep(NA,len2*len),ncol=len2)
  for(k in 1:len2)
  {
    for(j in 1:ll1)
      EUDn[,k][binSize==levs[j]]<-(v[binSize==levs[j],]%*%(doseBins[,j]^(1/n[k])))^n[k]
  }

maxEUD<-max(EUDn)
  list(EUD=EUDn/maxEUD
       ,maxEUD=maxEUD, n=n)
}

