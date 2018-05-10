llf<- function(beta,X,y,FUN) {
  #browser()
  sF2<-1e-17
  xbeta<- X%*%beta
  A<-ifelse(FUN(xbeta)==0, sF2,FUN(xbeta))
  A<-ifelse(FUN(xbeta)==1,1-sF2,A)
  neglog<--sum(y*log(A)+(1-y)*log(1-A))
  #if(is.nan(neglog) | is.na(neglog)|is.infinite(neglog))browser()
  neglog
}

#llf(start.Value,preDVH$EUD[,1],fractionation,toxicity,pnorm)
