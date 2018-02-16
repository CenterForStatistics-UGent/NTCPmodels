
llf<- function(param, EUD, frac,y,FUN) {
  xbeta<- ((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3])
  neglog<--sum(y*log(FUN(xbeta))+(1-y)*log(1-FUN(xbeta)))
  neglog
}
gra.lf<- function(param, EUD, frac,y,FUN) {
  #Not the most elegant way of doing things but it works
  if(names(formals(FUN)[2]) =="mean")FUN1<-dnorm
  if(names(formals(FUN)[2]) =="location")FUN1<-dlogis
  xbeta<- ((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3])
  d_m<--((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3]^2)
  d_TD50<--((((frac+param[1])/(2+param[1]))*(1/param[2]^2))*EUD)*(1/param[3])
  d_theta<-((((2-frac)/(2+param[1])^2)*(1/param[2]))*EUD)*(1/param[3])
  g_m<--sum(y*FUN(xbeta)/(FUN1(xbeta)*d_m)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_m))
  g_TD50<--sum(y*FUN(xbeta)/(FUN1(xbeta)*d_TD50)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_TD50))
  inter<-y*FUN(xbeta)/(FUN1(xbeta)*d_theta)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_theta)
  inter[frac==2]<-0
  g_theta<--sum(inter)
  grad<-c(g_theta,g_TD50,g_m)
  grad
}



NTCP<-function(DVH,
               fractionation,
               toxicity,
               link=c("probit","logit"),
               n=NULL)
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
    #d[,i][1:(lenn)]<-dvh[,1]
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
    n<-seq(0.01,1,length.out = len2)
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

  param <- c(theta = 3, TD50 = 20, m = 2)
  paras<-list()
  length(paras)<-len2
  link<-match.arg(link)
  switch (link,pnorm=probit,logit=plogis)

for(kk in 1:len2)
  {
  fit<- tryCatch(expr=optim(param, llf, gr=  gra.lf, EUD=EUDn[,kk],
                        frac=fractionation,y=toxicity,FUN=FUN,
                        method="L-BFGS-B", control=list(trace=TRUE, REPORT=1),
                        hessian=TRUE,
                        lower=c(0,0,-1000),upper=c(100,100,1000)),
             error=function(cond) {
               message(paste("EUD could not be computed for n =", n[kk]))
               return(NA)
             })
 paras[[kk]]<-fit

}
  out<-list(modelInfor=paras,dataInfo=EUDn,toxicity=toxicity,fractionation=fractionation,link=link)

}
fit<-NTCP(DVH,fractionation,toxicity,link = "logit",n=seq(0.01,1,length.out = 30)
          )
fit$link
xbeta<-function(param,EUD,fractionation)
{
(( (fractionation+param[1])/(2+param[1]))*(1/param[3])*(1/param[2])*EUD)-(1/param[1])
}

predict.NTCPmodels<-function(object, newDVH, n = NULL, type = c("link", "response",
                                                                 "class"))
{
  link<-xbeta(param=fit$modelInfor[[30]]$par,EUD=fit$dataInfo[,30], fractionation=fit$fractionation)
  if(fit$link=="logit")
  response<-plogis(link)
  class<-ifelse(response>0.5,0,1)

}



cv.NTCP<-function(DVH,
                  fractionation,
                  toxicity,
                  link=c("probit","logit"),
                  n=NULL,type.measure=c("AUC","Brier_scaled"),nfolds=5)