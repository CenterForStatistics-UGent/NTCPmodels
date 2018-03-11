auc.NTCPmodels<-function (y_pred, y_true)
{
  rank <- rank(y_pred)
  n_pos <- sum(y_true == 1)
  n_neg <- sum(y_true == 0)
  AUC <- (sum(rank[y_true == 1]) - n_pos * (n_pos + 1)/2)/(n_pos *
                                                             n_neg)
  return(AUC)
}

cv.NTCP<-function(DVH,
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
    callNTCP<-NTCP(DVH[-testID],
                   fractionation[-testID],
                   toxicity[-testID],
                   link,n)
    pred<-predict.NTCPmodels(callNTCP,DVH[testID],fractionation[testID],type="response")
    measurefolds[[i]]<-apply(pred,2,function(x)auc.NTCPmodels(x,toxicity[testID]))
    }
DT<-data.frame(unlist(measurefolds))
out<-aggregate(DT,list(names(unlist(measurefolds))),"mean")
names(out)<-c("n",measure)
class(out)<-"cvNTCPmodels"
out
    }


