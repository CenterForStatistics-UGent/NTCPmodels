#' Fits NTCP models.
#'
#' ```NTCP```
#'
#' This function fits NTCP, Lyman models, and can account for fractionation. Logit  models
#' can also be fitted
#'
#'
#'
#' @param DVH  A list of DVH data for all patients. Each element in the list is  a two columns ```data.frame``` of DVH
#' data for a patient. The first column of the the ```data.frame```` must be the dose and the second
#' column must be the  volume recieving atmost that dose.
#' @param fractionation A vector of same length as the ```DVH``` list,  which contains the dose per
#' fraction for each patient.
#' @param XClin A numeric matrix of clinical variables
#' with each row corresponding to a patient and the columns are the variables.
#' @param  toxicity A vector of same length as the fractionation, contain toxicity status of
#' all patients. Toxicity should be coded as 1 and no-toxicity coded as 0
#' @param link The type of link function,should either be "logit" or "probit" link.
#' @param  n A vector or numeric constant, representing the ```n``` parameter in the
#' Kutcher and Burman EUD reduction scheme, and should lie between 0 and 1.
#' The defualt value is ```NULL```. If no ```n``` constant or vector is supplied, NTCP will
#' try a grid of ```n``` values between 0 and 1 and return the fitted model parameters for each
#' ```n```.
#' @return
#' @examples
#' fit<-NTCP(DVH=DVH,fractionation=fractionation, toxicity=toxicity,link="logit")
#' summary(fit)



NTCP<-function(DVH,
               XClin=NULL,
               fractionation,
               toxicity,
               link=c("probit","logit"),
               n=NULL)
{
  link=match.arg(link)
  if(!link %in% c("probit","logit"))
    stop("Link must be either 'probit' or  'logit'")
    pClin<-dim(XClin)[2]

  preDVH<-prepare.DVH(DVH,n=n)
  lln<-length(preDVH$EUD[,1])
  n<-preDVH$n
  len2<-length(preDVH$n)
  if(!is.null(XClin))
  {

    names(XClin)<-paste("X_",c(1:pClin),sep = "")
    paras<-data.frame(matrix(rep(NA,len2*(4+pClin)),ncol=(4+pClin)))
    colnames(paras)<-c("n", "beta0","beta1","beta3",names(XClin))
  }
  else
  {
    paras<-data.frame(matrix(rep(NA,len2*(4)),ncol=4))
    colnames(paras)<-c("n", "beta0","beta1","beta3")
  }

  paras$n<-round(n,2)

  if(missing(fractionation))
  fractionation<-rep(2,length(toxicity))

  if(is.null(XClin))
  {
    beta_0 =c(0,0,0)
    low<-  c(-100,0,0)
    up<-c(-0.0001,100,100)
  }
  else
  {
    beta_0 =c(0,0,0,rep(0,dim(XClin)[2]))
    low<-c(-100,0,0,rep(-100,pClin))
    up<-c(-0.0001,100,100,rep(100,pClin))
  }


for(kk in 1:len2)
  {
  X<-cbind(rep(1,lln),preDVH$EUD[,kk],fractionation*preDVH$EUD[,kk],XClin)
  X<-as.matrix(X)
  initialMod<- stats::optim(fn= function(beta,X,y){
                                  sum((y -X%*%beta)^2)
                                     },
                     lower = low,
                     upper=up,
                     y=toxicity,X=X,
                     par = beta_0,
                     method="L-BFGS-B")


  start.Value<-initialMod$par
  if(link=="probit")
  fit<- tryCatch(expr=stats::optim(start.Value, llf,X=X,
                            #gr=  gra.lf,
                        y=toxicity,FUN=stats::pnorm,
                        method="L-BFGS-B", control=list(trace=TRUE, REPORT=1),
                        hessian=TRUE,
                        lower=low,upper=up
                        ),
             error=function(cond) {
               message(paste("EUD could not be computed for n =", n[kk]))
               return(NA)
             })
  else
    fit<- tryCatch(expr=stats::optim(start.Value, llf,X=X,
                              #gr=  gra.lf,
                              y=toxicity,FUN=stats::plogis,
                              method="L-BFGS-B",
                              control=list(trace=TRUE, REPORT=1),
                              hessian=TRUE,
                              lower=low,upper=up
                              ),
                   error=function(cond) {
                     message(paste("EUD could not be computed for n =", n[kk]))
                     return(NA)
                   })
  if(!is.na(fit[[1]][1]))
    paras[kk,-1]<-fit$par
  else
    paras[kk,-1]<-rep(NA,dim(X)[2])


}
paras<-stats::na.omit(paras)

  # exclude<-apply(paras[,c(2:4)],1,function(x)
  # {
  #
  #   ret<-1
  #   if(sum(x==start.Value)==3)
  #   ret<-0
  #   ret
  # }
  # )

  m<--as.numeric(1/paras[,2])
  theta<-ifelse(paras[,4]==0,100,as.numeric(paras[,3]/paras[,4]))
  if(sum(paras[,4]==0)<=sum(paras[,3]==0))
    TD_50<-ifelse(paras[,4]==0,100,as.numeric(preDVH$maxEUD/(m*(2+theta)*paras[,4])))
  if(sum(paras[,4]==0)>sum(paras[,3]==0))
    TD_50<-ifelse(paras[,3]==0,100,as.numeric((theta*preDVH$maxEUD)/(m*(2+theta)*paras[,3])))


  if(dim(paras)[2]==4)
  paras1<-cbind(paras[,1],m, TD_50,theta)
  else
    paras1<-cbind(paras[,1],m, TD_50,theta,paras[,-c(1:4)])
  EUDout<- as.data.frame(preDVH$EUD)
  names(EUDout)<-paste("EUD_",sep="",round(n,2))
  out<-list(lymanPara=paras1,
            betas=paras[,-1],
            dataInfo=EUDout,
            XClinical=XClin,
            toxicity=toxicity,
            fractionation=fractionation,
            link=link,
            n=n)
  class(out)<-"NTCPmodels"
  out

}

#Looklihood function

#Gradient function


