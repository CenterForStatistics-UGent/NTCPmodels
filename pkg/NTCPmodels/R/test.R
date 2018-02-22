fit<-NTCP(DVH,fractionation,toxicity)
cv.fitUI<-cv.NTCP(DVH,fractionation,toxicity)
cv.fitUF<-cv.NTCP(DVH,fractionation,toxicity,link = "logit")



par(mfrow=c(1,2))
plot.cvNTCPmodels(cv.fitUI)
fitUI<-NTCP(DVH,fractionation,toxicity)
#n= 0.2  theta=29.719783  TD_50=48.75401  m=38.449429

plot.cvNTCPmodels(cv.fitUF)
fitUF<-NTCP(DVH,fractionation,toxicity,main="Frequency")
#n=0.39 theta=28.001573  TD_50=47.50115  m=29.680333