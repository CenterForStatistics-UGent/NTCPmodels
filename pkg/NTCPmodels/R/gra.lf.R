# gra.lf<- function(param, Xy,FUN) {
#     sF<-1e5
#     sF2<-1e-17
#     X[,-1]<-X[,-1]/sF
#     grad <- beta*0
#     xbeta<-X%*%beta
#     A<-ifelse(FUN(xbeta)==0,sF2,FUN(xbeta))
#     A<-ifelse(FUN(xbeta)==1,1-sF2,A)
#     if(names(formals(FUN)[2]) =="mean")pdf<-stats::dnorm
#     if(names(formals(FUN)[2]) =="location")pdf<-stats::dlogis
#     for (k in 1:K) {
#       grad[k] <- sum(y*(pdf(xbeta)*X[,k])/(A)+ (y-1)*(pdf(xbeta)*X[,k])/(1-A))
#     }
#     grad[-1]<-grad[-1]*sF
#     return(-grad)
# }
