#' Fits NTCP models.
#'
#' \code{NTCP}
#'
#' This function fits NTCP, Lyman models, and can account for fractionation. Logit  models
#' can also be fitted
#'
#'
#'
#' @param DVH  A list of DVH data for all patients. Each element in the list is  a two columns \code{data.frame} of DVH
#' data for a patient. The first column of the the \code{data.frame} must be the dose and the second
#' column must be the  volume recieving atmost that dose.
#' @param fractionation A vector of same length as the \code{DVH} list,  which contains the dose per
#' fraction for each patient.
#' @param  toxicity A vector of same length as the fractionation, contain toxicity status of
#' all patients. Toxicity should be coded as 1 and no-toxicity coded as 0
#' @param  n A vector or numeric constant, representing the \code{n} parameter in the
#' Kutcher and Burman EUD reduction scheme, and should lie between 0 and 1.
#' The defualt value is \code{NULL}. If no \code{n} constant or vector is supplied, NTCP will
#' try a grid of \code{n} values between 0 and 1 and return the fitted model parameters for each
#' \code{n}.
#' @param start.Value The starting values used by the \code{optim} function fit
#' the model paramters (\code{theta},\code{TD_50},\code{m}). Must be a vector of length three,
#' the first element must be \code{theta}, second \code{TD_50} and third \code{m}.
#' The default starting values are \code{c(theta = 4, TD50 = 30, m = 1)}
#' @return
#' @examples



NTCP<-function(DVH,
               fractionation,
               toxicity,
               link=c("probit","logit"),
               n=NULL,
               start.Value=c(theta = 4, TD50 = 30, m = 1))
{
  link=match.arg(link)
  if(!link %in% c("probit","logit"))
    stop("Link must be either 'probit' or  'logit'")
  preDVH<-prepare.DVH(DVH,n=n)

  n<-preDVH$n
  len2<-length(preDVH$n)
  paras<-data.frame(matrix(rep(NA,len2*4),ncol=4))
  colnames(paras)<-c("n", "theta","TD_50","m")
  paras$n<-round(n,2)


for(kk in 1:len2)
  {
  if(link=="probit")
  fit<- tryCatch(expr=optim(start.Value, llf, gr=  gra.lf, EUD=preDVH$EUD[,kk],
                        frac=fractionation,y=toxicity,FUN=pnorm,
                        method="L-BFGS-B", control=list(trace=TRUE, REPORT=1),
                        hessian=TRUE,
                        lower=c(0,0,-1000),upper=c(100,100,1000)),
             error=function(cond) {
               message(paste("EUD could not be computed for n =", n[kk]))
               return(NA)
             })
  else
    fit<- tryCatch(expr=optim(start.Value, llf, gr=  gra.lf, EUD=preDVH$EUD[,kk],
                              frac=fractionation,y=toxicity,FUN=plogis,
                              method="L-BFGS-B", control=list(trace=TRUE, REPORT=1),
                              hessian=TRUE,
                              lower=c(0,0,-1000),upper=c(100,100,1000)),
                   error=function(cond) {
                     message(paste("EUD could not be computed for n =", n[kk]))
                     return(NA)
                   })
  if(!is.na(fit[[1]][1]))
    paras[kk,2:4]<-fit$par
  else
    paras[kk,2:4]<-rep(NA,3)


}
  out<-list(modelInfor=paras,dataInfo=preDVH$EUD,toxicity=toxicity,fractionation=fractionation,link=link,n=n)

}

#Looklihood function
llf<- function(param, EUD, frac,y,FUN) {
  xbeta<- ((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3])
  neglog<--sum(y*log(FUN(xbeta))+(1-y)*log(1-FUN(xbeta)))
  neglog
}
#Gradient function
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

  list(EUD=EUDn,n=n)
}
