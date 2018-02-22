llf<- function(param, EUD, frac,y,FUN) {
  xbeta<- ((((frac+param[1])/(2+param[1]))*(1/param[2]))*EUD-1)*(1/param[3])
  neglog<--sum(y*log(FUN(xbeta))+(1-y)*log(1-FUN(xbeta)))
  neglog
}

#llf(start.Value,preDVH$EUD[,1],fractionation,toxicity,pnorm)
