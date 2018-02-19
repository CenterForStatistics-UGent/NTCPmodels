library(glmnet)
xbeta<-function(param,EUD,fractionation)
{
 param<-as.numeric(param)
  (( (fractionation+param[1])/(2+param[1]))*(1/param[3])*(1/param[2])*EUD)-(1/param[3])
}

predict.NTCPmodels<-function (object, newDVH,newFractionation, type = c("link", "response","class"))
{
  type = match.arg(type)
  len1<-length(object$n)
  nNotNA<-c(1:len1)[!is.na(object$modelInfor[,2])]

  len3<-length(nNotNA)


  if (missing(newDVH))
  {
    EUD<-object$dataInfo
    len2<-dim(EUD)[1]
    Link<-matrix(rep(NA,len2*len3),ncol = len3)
    for(i in 1:len3)
    {
      Link[,i]<-xbeta(object$modelInfor[nNotNA[i],2:4],
                      EUD[,nNotNA[i]],
                      object$fractionation)
    }
    pred<-Link
  }
  else
  {
    if(missing(newFractionation))
      stop("Missing newFractionation data")
    if(length(newDVH)!=length(newFractionation))
      stop("difference in length  of newDVH and length of newFractiontion")
    prepDVH<-prepare.DVH(newDVH,n=object$n[nNotNA])
    EUD<-prepDVH$EUD
    len2<-dim(EUD)[1]
    Link<-matrix(rep(NA,len2*len3),ncol = len3)
    for(i in 1:len3)
    {
      Link[,i]<-xbeta(object$modelInfor[nNotNA[i],2:4],
                      EUD[,i],
                      newFractionation)
    }
    pred<-Link
  }



  if(type %in% c("response","class"))
  {
    if(object$link=="probit")
      response<-apply(Link,2,pnorm)
    else
      response<-apply(Link,2,plogis)
    pred<-response
  }


  if(type=="class")
  {
    class<-ifelse(response<0.5,0,1)
    pred<-class
  }

  pred<-data.frame(pred)
  names(pred)<-round(object$n[nNotNA],2)
  pred
}


