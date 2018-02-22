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

if(missing(fractionation))
  fractionation<-rep(2,length(toxicity))

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
                              lower=c(0,0,-1000),upper=c(100,80,1000)),
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
  class(out)<-"NTCPmodels"
  out

}

#Looklihood function

#Gradient function


