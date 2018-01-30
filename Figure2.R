library(deSolve)
library(minqa)

GrowthTime <- 36
NAD.c <- 168.62
NADP.c <- 2.442

NAD_Measurement <- c(0.246821255,0.379819763,0.470399904,0.610853991,0.709800956,0.863714214,0.944326188)
NADP_Measurement <- c(0.177051892,0.260108406,0.402180048,0.56771244,0.677802989,0.810994031,0.909691225)

NAD_NADP_Measurement <- cbind(NAD_Measurement,NADP_Measurement)

NAD_ODE <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dNAD <- (1-NAD)*F1/NAD.c
    dNADP <- (NAD-NADP)*F2/NADP.c
    return(list(c(dNAD,dNADP)))
  })
}

yini <- c(NAD = 0, NADP = 0)
times <- seq(0, GrowthTime, by = 0.1)

NAD <- function(f) {
  pars <- c(F1=f[1], F2=f[2])
  out <- ode(func = NAD_ODE, y = yini, parms = pars, times = times,method="lsoda")
  return(out[c(41,61,81,116,156,241,361),2:3])
}

Res <- function(f) {
  sum((NAD_NADP_Measurement-NAD(f))^2)
}

system.time(fluxes<-bobyqa(c(0.02,0.01),lower=c(0,0),fn=Res))