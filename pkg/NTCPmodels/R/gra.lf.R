gra.lf<- function(param, EUD, frac,y,FUN) {
  #Not the most elegant way of doing things but it works
  if(names(formals(FUN)[2]) =="mean")FUN1<-dnorm
  if(names(formals(FUN)[2]) =="location")FUN1<-dlogis
  xbeta<- ((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3])
  d_m<--((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3]^2)
  d_TD50<--((((frac+param[1])/(2+param[1]))*(1/param[2]^2))*EUD)*(1/param[3])
  d_theta<-((((2-frac)/(2+param[1])^2)*(1/param[2]))*EUD)*(1/param[3])
  g_m<--sum(y*FUN(xbeta)/(FUN1(xbeta)*d_m)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_m))
  g_TD50<--sum(y*FUN(xbeta)/(FUN1(xbeta)*d_TD50)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_TD50))
  inter<-y*FUN(xbeta)/(FUN1(xbeta)*d_theta)-(1-y)*(1-FUN(xbeta))/(FUN1(xbeta)*d_theta)
  inter[frac==2]<-0
  g_theta<--sum(inter)
  grad<-c(g_theta,g_TD50,g_m)
  grad
}
