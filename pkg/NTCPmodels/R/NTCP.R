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
#' @param link The type of link function,should either be "logit" or "probit" link.
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
#' fit<-NTCP(DVH,fractionation, toxicity,link="logit")
#' summary(fit)



NTCP<-function(DVH,
               fractionation,
               toxicity,
               link=c("probit","logit"),
               n=NULL)
{
  link=match.arg(link)
  if(!link %in% c("probit","logit"))
    stop("Link must be either 'probit' or  'logit'")
  preDVH<-prepare.DVH(DVH,n=n)
  lln<-length(preDVH$EUD[,1])
  n<-preDVH$n
  len2<-length(preDVH$n)
  paras<-data.frame(matrix(rep(NA,len2*4),ncol=4))
  colnames(paras)<-c("n", "beta0","beta1","beta3")
  paras$n<-round(n,2)
if(missing(fractionation))
  fractionation<-rep(2,length(toxicity))

low<-  c(-100,0,0)
up<-c(0,100,100)
for(kk in 1:len2)
  {
  X<-cbind(rep(1,lln),preDVH$EUD[,kk],fractionation*preDVH$EUD[,kk])
  beta_0 =c(0,0,0)
  initialMod<- optim(fn= function(beta,X,y){
                                  sum((y -X%*%beta)^2)
                                     },
                     lower = low,
                     upper=up,
                     y=toxicity,X=X,
                     par = beta_0,
                     method="L-BFGS-B")


  start.Value<-initialMod$par
  if(link=="probit")
  fit<- tryCatch(expr=optim(start.Value, llf,X=X,
                            #gr=  gra.lf,
                        y=toxicity,FUN=pnorm,
                        method="L-BFGS-B", control=list(trace=TRUE, REPORT=1),
                        hessian=TRUE,
                        lower=low,upper=up
                        ),
             error=function(cond) {
               message(paste("EUD could not be computed for n =", n[kk]))
               return(NA)
             })
  else
    fit<- tryCatch(expr=optim(start.Value, llf,X=X,
                              #gr=  gra.lf,
                              y=toxicity,FUN=plogis,
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
    paras[kk,2:4]<-fit$par
  else
    paras[kk,2:4]<-rep(NA,3)


}


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
  theta<-as.numeric(paras[,3]/paras[,4])
  TD_50<-as.numeric(preDVH$maxEUD/(m*(2+theta)*paras[,4]))
  paras<-cbind(paras[,1],m, TD_50,theta)
  EUDout<- as.data.frame(preDVH$EUD)
  names(EUDout)<-paste("EUD_",sep="",round(n,2))
  out<-list(modelInfor=paras,dataInfo=EUDout,toxicity=toxicity,fractionation=fractionation,link=link,n=n)
  class(out)<-"NTCPmodels"
  out

}

#Looklihood function

#Gradient function


