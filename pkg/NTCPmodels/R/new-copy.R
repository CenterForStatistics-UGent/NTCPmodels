



computef <- function(var, param, J = J,EQD2=TRUE) {
  theta<-0
  if(EQD2)
  {
    if(length(param)!=4)stop("parameters are set wrongly")
    theta<-param[1]
    n<-param[2]
    D50<-param[3]
    m<-param[4]
  }
  else
  {
    if(length(param)!=3)stop("parameters are set wrongly")
    n<-param[1]
    D50<-param[2]
    m<-param[3]
  }
  v <- var[((2 * J) + 1):((2 * J) + J)]
  D <- var[1:J]

  d<-2
  if(EQD2)
  d <- var[(J + 1):(2 * J)]

  summand1<- (((D * d) + (D * theta))/(2 +theta))
  if(EQD2)
  summand11<-(((2*D) -(D * d))/((2 +theta)^2))

  sum1 <- as.numeric(v) %*%matrix(summand1^(1/n))
  sum11 <-as.numeric(v)%*%matrix(as.numeric((summand1^(1/n))*log(summand1)))
  if(EQD2)
  sum111<-as.numeric(v)%*%matrix((1/n)*(summand1^((1-n)/n))*as.numeric(summand11))

  f<-(1/m)*((sum1^n)*(1/D50) - 1)

  f_dm<--(1/m^2)*((sum1^n)*(1/D50) - 1)
  f_D50<--(1/m)*((sum1^n)*(1/D50^2))
  f_n<-(1/m)*(1/D50)*(sum1^n)*log(sum1)*sum11

  if(EQD2)
  f_theta<-(1/m)*(1/D50)*n*(sum1^(n-1))*sum111

  if(EQD2)
  return(list(f=f,f_dm=f_dm,f_D50=f_D50,f_n=f_n,f_theta=f_theta))
  if(!EQD2)
  return(list(f=f,f_dm=f_dm,f_D50=f_D50,f_n=f_n))
}
logit.lf<- function(param, datafit, J, y,EQD2) {
  fit <- apply(datafit, 1, function(x,param,J,EQD2){
    computef(x, param, J = J,EQD2)$f
  },param = param,J = J,EQD2)
  prob1 <- plogis(fit)
  logexb <- log(prob1)
  y0 <- 1 - y
  logexb0 <- log(1 - prob1)
  yt <- t(y)
  y0t <- t(y0)
  logl <- -sum(yt%*%logexb + y0t%*%logexb0)
  return(logl)
}
logit.gr <-function(param, datafit, J,y,EQD2) {
  fit <- apply(datafit, 1, function(x,param,J,EQD2){
    unlist(computef(x, param, J = J,EQD2))
  },param = param,J = J,EQD2)

  fit<-data.frame(t(fit))
  exb <- exp(fit$f)
  prob <- exb/(1+exb)
  g_m<-sum(y*fit$f_dm-prob)
  g_D50<-sum(y*fit$f_D50-prob)
  g_n<-sum(y*fit$f_n-prob)

  if(EQD2)
  {
    g_theta<-sum(y*fit$f_theta-prob)
    grad<-c(g_theta,g_n,g_D50,g_m)
  }
  if(!EQD2)
  {
    grad<-c(g_n,g_D50,g_m)
  }
  return(-grad)

}
#
# param <- c(theta =2, n = 2, D50 = 60, m = 2)
# datafit <- ntcpObject$datafit
# J <- ntcpObject$J
# y <- ntcpObject$tox
#
# param <- c(theta=2,n =1, D50 = 20, m = 12)
# fit <- apply(datafit, 1, function(x,param,J,EQD2){
#  unlist(computef(x, param, J = J,EQD2))
# },param = param,J = J,EQD2=TRUE)
#
# t(fit)
# ff<-data.frame(t(fit))$f
# yy<-ifelse(plogis(ff)>0.51,1,0)
# table(yy)
# param <- c(n = 1, D50 = 60, m = 1)
# logit.lf(param, datafit, J,yy,EQD2 = FALSE)
# logit.gr(param, datafit, J,yy,EQD2 = FALSE)
#
#
#
# #
#
# # # logit gradient function ##
# # param <- c(n = 1.1, D50 = 50, m = 3)
# # logitmodel <- optim(param, logit.lf, gr=logit.gr, datafit=datafit,
# #                     J=J,y=yy,EQD2=FALSE,
# #                     method="BFGS", control=list(trace=TRUE, REPORT=1),
# #                     hessian=TRUE)
# # coeffs <- logitmodel$par
# # covmat <- solve(logitmodel$hessian)
# # stderr <- sqrt(diag(covmat))
# # zscore <- coeffs/stderr
# # pvalue <- 2*(1 - pnorm(abs(zscore)))
# # results <- cbind(coeffs,stderr,zscore,pvalue)
# # colnames(results) <- c("Coeff.", "Std. Err.", "z", "p value")
# # print(results)
# # #
# #
# # ## Graph likelihood for age coefficient, holding all else at maximum
# #
