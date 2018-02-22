xbeta<-function(param,EUD,fractionation)
{
  param<-as.numeric(param)
  (( (fractionation+param[1])/(2+param[1]))*(1/param[3])*(1/param[2])*EUD)-(1/param[3])
}
