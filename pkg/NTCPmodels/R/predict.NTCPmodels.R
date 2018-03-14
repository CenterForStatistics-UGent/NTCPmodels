
predict.NTCPmodels<-function (object, newDVH,newXClin=NULL,newFractionation, type = c("link", "response","class"))
{
  type = match.arg(type)
  len1<-length(object$n)
  nNotNA<-c(1:len1)[!is.na(object$betas[,2])]

  len3<-length(nNotNA)

  #Using object info
  if (missing(newDVH))
  {
    EUD<-object$dataInfo
    len2<-dim(EUD)[1]
    Link<-matrix(rep(NA,len2*len3),ncol = len3)
    for(i in 1:len3)
    {
      X<-cbind(rep(1,len2),object$dataInfo[,i],object$dataInfo[,i]*object$fractionation,object$XClinical)
      Link[,i]<-as.matrix(X)%*%as.numeric(object$betas[i,])
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
      X<-cbind(rep(1,len2), EUD[,i], EUD[,i]*newFractionation,newXClin)
      Link[,i]<-as.matrix(X)%*%as.numeric(object$betas[i,])
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


