auc.NTCPmodels<-function (y_pred, y_true)
{
  rank <- rank(y_pred)
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  AUC <- (sum(rank[y_true == 1]) - n_pos * (n_pos + 1)/2)/(n_pos *
                                                             n_neg)
  return(AUC)
}

cv.NTCP<-function(DVH,XClin=NULL,
          fractionation,
          toxicity,n= NULL,
          link=c("probit","logit"),
          measure = "auc",
          nfolds = 10, foldid)
{
  #type.measure<-match.arg(type.measure)

  N<-length(DVH)
  folds <- cut(seq(1,N),breaks=nfolds,labels=FALSE)
  folds<- folds[sample(1:N)]
  measurefolds<-list()
  length(measurefolds)<-nfolds
  for(i in 1:nfolds){

    testID <- which(folds==i,arr.ind=TRUE)
    callNTCP<-NTCP(DVH=DVH[-testID],
                   XClin=XClin[-testID,],
                   fractionation=fractionation[-testID],
                   toxicity=toxicity[-testID],
                   link,n)
    pred<-predict.NTCPmodels(callNTCP,newDVH=DVH[testID],newXClin=XClin[testID,],newFractionation=fractionation[testID],type="response")
    measurefolds[[i]]<-apply(pred,2,function(x)auc.NTCPmodels(x,toxicity[testID]))
  }
  tmp<-unlist(measurefolds)
  rmv<-!is.nan(tmp)
  tmp<-tmp[rmv]
DT<-data.frame(tmp)
out<-aggregate(DT,list(names(tmp)),"mean")
names(out)<-c("n","measure")
if(measure=="auc")
{
  AUC<-out[,2]
  N1<-table(toxicity)[1]
  N2<-table(toxicity)[2]

  Q1<-AUC/(2-AUC)
  Q2<-(2*AUC^2)/(1+AUC)

  seAUC<-sqrt((AUC*(1-AUC)+(N1-1)*(Q1-AUC^2)+(N2-1)*(Q2-AUC^2))/(N1*N2))
  sd<-seAUC
}
out<-list(out,sd=sd,type.measure=measure)
class(out)<-"cvNTCPmodels"
out
    }


