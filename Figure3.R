library(deSolve)
library(minqa)

GrowthTime <- 36

NAD_Measurement <- c(0.246821255,0.379819763,0.470399904,0.610853991,0.709800956,0.863714214,0.944326188)
NADP_Measurement <- c(0.177051892,0.260108406,0.402180048,0.56771244,0.677802989,0.810994031,0.909691225)

NAD_NADP_Measurement <- cbind(NAD_Measurement,NADP_Measurement)

NAD_ODE <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dNAD <- (1-NAD)*F1
    dNADP <- (NAD-NADP)*F2
    return(list(c(dNAD,dNADP)))
  })
}
yini <- c(NAD = 0, NADP = 0)
times <- seq(0, GrowthTime, by = 0.1)

NAD <- function(f) {
  #  pars <- c(PCV, Seeding, FGlc=f[1], FGlySer=f[2], FConSer=(SerNetConsumption+f[1]+f[2]),FSerGly=f[3],FConGly=GlyNetConsumption+f[3])
  pars <- c(F1=f[1], F2=f[2])
  out <- ode(func = NAD_ODE, y = yini, parms = pars, times = times,method="lsoda")
  return(out[c(41,61,81,116,156,241,361),2:3])
}

Res <- function(f) {
  sum((NAD_NADP_Measurement-NAD(f))^2)
}

#system.time(print(optim(c(10,10,10,10,10),fn=Res,control=list(maxit=10000))))
system.time(fluxes<-bobyqa(c(0.02,0.01),lower=c(0,0),fn=Res))

#Test <- expand.grid(In = seq(6, 13, len = 701), Out = seq(5, 5, len = 1))

#for (i in 1:nrow(Test)) {
#  Test[i,3] <- sum((NAD_Measurement-NAD(c(Test[i,1],Test[i,2])))^2)
#}



#FGlc=5.83
#FGlySer=15.74
#FConSer=68.51
#FSerGly=40.25
#FConGly=41.14

times <- seq(0, GrowthTime, by = 0.01)
#pars <- c(PCV, Seeding,FGlc=15.927,FGlySer=27.333,FConSer=82.314,FSerGly=38.107,FConGly=46.580)
pars <- c(F1=fluxes$par[1], F2=fluxes$par[2])
system.time(output <- ode(func = NAD_ODE, y = yini, parms = pars, times = times))
outputM <- data.frame(output)

write.csv(outputM,file="NAD_NADP.csv")

outputM$NAD <- (outputM$NAD_U+outputM$NAD_L)/outputM$Cell
outputM$NAD_UF <- outputM$NAD_U/(outputM$NAD_U+outputM$NAD_L)
outputM$NAD_LF <- outputM$NAD_L/(outputM$NAD_U+outputM$NAD_L)

matplot(output[,"time"], outputM[,6:7], type = "l", xlab = "time", ylab = "Concentration(uM)",main = "NAD", lwd = 2)

