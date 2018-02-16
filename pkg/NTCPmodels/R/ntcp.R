#'Fits NTCP models from an ntcpObject object
#'
#'ntcp returns a model object with fitted parameter values.
#'This function takes in an ntcpObject object, and can fit four different NTCP models.
#'@param ntcpObject an ntcp object created by the function createNtcpObject.
#'@param model the type model to fit it should be one of: LymanEUD, LymanEUDEQD2, LogitEUD
#'and LogitEUDEQD2. The default is LymanEUD.
#'@param nstart the number of starting values from which the fitting engine will use as
#'the starting. If the defualt is nstart=5.
#'@examples
#'@export
ntcp<-function(ntcpObject,model=c("LymanKB.EUD",
                                  "LymanKB.EUD.EQD2",
                                  "LogitKB.EUD",
                                  "LogitKB.EUD.EQD2"),
                                 nstart=5){
  #check to makge sure ntcpObject is an ntcp object
  datafit<-ntcpObject$datafit
  J<-ntcpObject$J
  tox<-ntcpObject$tox
  model<-match.arg(model)
  switch (model,
          LymanKB.EUD=fitEUD(datafit,tox=tox,N=nstart,link = "probit",J=J),
          LymanKB.EUD.EQD2=fitEUDEQD2(datafit,tox=tox,N=nstart,link = "probit",J=J),
          LogitKB.EUD=fitEUD(datafit,N=nstart,tox=tox,link = "logit",J=J),
          LogitKB.EUD.EQD2=fitEUDEQD2(datafit,tox=tox,N=nstart,link = "logit",J=J)

  )
}





#fit specific models
