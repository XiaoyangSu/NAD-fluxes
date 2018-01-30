require(deSolve)
require(minqa)
require(numDeriv)
require(MASS)

NADNAMmod0D1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
####NAM tracer####
    dNAM_S0.NAM <- 3.5066*Time-1.7139
    dNAM_S3.NAM <- -0.001*Time+0.0245
    dNAM_S4.NAM <- -3.5066*Time+1.7139+0.001*Time-0.0245
    dNA_0.NAM <- 0
    dNA_3.NAM <- 0
    dNA_4.NAM <- 0
    dTrp_0.NAM <- 0
    dTrp_3.NAM <- 0
    dTrp_4.NAM <- 0
    dNAD_0.NAM <- (f1*(Trp_0.NAM-NAD_0.NAM)+f2*(NA_0.NAM-NAD_0.NAM)+
                     f3*(NAM_0.NAM-NAD_0.NAM))/cNAD
    dNAD_3.NAM <- (f1*(Trp_3.NAM+Trp_4.NAM-NAD_3.NAM)+
                     f2*(NA_3.NAM+NA_4.NAM-NAD_3.NAM)+
                     f3*(NAM_3.NAM+NAM_4.NAM-NAD_3.NAM))/cNAD
    dNAM_0.NAM <- ((f1+f2+f3)*(NAD_0.NAM-NAM_0.NAM)+f4*(NAM_S0.NAM-NAM_0.NAM))/cNAM
    dNAM_3.NAM <- ((f1+f2+f3)*(NAD_3.NAM-NAM_3.NAM)+f4*(NAM_S3.NAM-NAM_3.NAM))/cNAM
    dNAM_4.NAM <- ((f1+f2+f3)*(0-NAM_4.NAM)+f4*(NAM_S4.NAM-NAM_4.NAM))/cNAM
####NA tracer####    
    dNAM_S0.NA <- -0.002
    dNAM_S6.NA <- 0.002
    dNA_0.NA <- 0
    dNA_6.NA <- 0
    dTrp_0.NA <- 0
    dTrp_6.NA <- 0
    dNAD_0.NA <- (f1*(Trp_0.NA-NAD_0.NA)+f2*(NA_0.NA-NAD_0.NA)+
                    f3*(NAM_0.NA-NAD_0.NA))/cNAD
    dNAD_6.NA <- (f1*(Trp_6.NA-NAD_6.NA)+
                     f2*(NA_6.NA-NAD_6.NA)+
                     f3*(NAM_6.NA-NAD_6.NA))/cNAD
    dNAM_0.NA <- ((f1+f2+f3)*(NAD_0.NA-NAM_0.NA)+f4*(NAM_S0.NA-NAM_0.NA))/cNAM
    dNAM_6.NA <- ((f1+f2+f3)*(NAD_6.NA-NAM_6.NA)+f4*(NAM_S6.NA-NAM_6.NA))/cNAM
####Tryptophan tracer#### 
    dNAM_S0.Trp <- -0.005*Time-0.0043
    dNAM_S6.Trp <- 0.005*Time+0.0043
    dTrp_0.Trp <- 0
    dTrp_6.Trp <- 0
    dNA_0.Trp <- 0
    dNA_6.Trp <- 0
    dNAD_0.Trp <- (f1*(Trp_0.Trp-NAD_0.Trp)+f2*(NA_0.Trp-NAD_0.Trp)+
                     f3*(NAM_0.Trp-NAD_0.Trp))/cNAD
    dNAD_6.Trp <- (f1*(Trp_6.Trp-NAD_6.Trp)+f2*(NA_6.Trp-NAD_6.Trp)+
                     f3*(NAM_6.Trp-NAD_6.Trp))/cNAD
    dNAM_0.Trp <- ((f1+f2+f3)*(NAD_0.Trp-NAM_0.Trp)+f4*(NAM_S0.Trp-NAM_0.Trp))/cNAM
    dNAM_6.Trp <- ((f1+f2+f3)*(NAD_6.Trp-NAM_6.Trp)+f4*(NAM_S6.Trp-NAM_6.Trp))/cNAM
    return(list(c(dNAM_S0.NAM,dNAM_S3.NAM,dNAM_S4.NAM,
                  dNA_0.NAM,dNA_3.NAM,dNA_4.NAM,
                  dTrp_0.NAM,dTrp_3.NAM,dTrp_4.NAM,
                  dNAD_0.NAM,dNAD_3.NAM,dNAM_0.NAM,dNAM_3.NAM,dNAM_4.NAM,
                  dNAM_S0.NA,dNAM_S6.NA,dNA_0.NA,dNA_6.NA,dTrp_0.NA,dTrp_6.NA,
                  dNAD_0.NA,dNAD_6.NA,dNAM_0.NA,dNAM_6.NA,
                  dNAM_S0.Trp,dNAM_S6.Trp,dNA_0.Trp,dNA_6.Trp,dTrp_0.Trp,dTrp_6.Trp,
                  dNAD_0.Trp,dNAD_6.Trp,dNAM_0.Trp,dNAM_6.Trp)))
  })
}

NADNAMmod0D2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ####NAM tracer####
    dNAM_S0.NAM <- -0.059/Time
    dNAM_S3.NAM <- -0.001*Time+0.0245
    dNAM_S4.NAM <- 0.059/Time+0.001*Time-0.0245
    dNA_0.NAM <- 0
    dNA_3.NAM <- 0
    dNA_4.NAM <- 0
    dTrp_0.NAM <- 0
    dTrp_3.NAM <- 0
    dTrp_4.NAM <- 0
    dNAD_0.NAM <- (f1*(Trp_0.NAM-NAD_0.NAM)+f2*(NA_0.NAM-NAD_0.NAM)+
                     f3*(NAM_0.NAM-NAD_0.NAM))/cNAD
    dNAD_3.NAM <- (f1*(Trp_3.NAM+Trp_4.NAM-NAD_3.NAM)+
                     f2*(NA_3.NAM+NA_4.NAM-NAD_3.NAM)+
                     f3*(NAM_3.NAM+NAM_4.NAM-NAD_3.NAM))/cNAD
    dNAM_0.NAM <- ((f1+f2+f3)*(NAD_0.NAM-NAM_0.NAM)+f4*(NAM_S0.NAM-NAM_0.NAM))/cNAM
    dNAM_3.NAM <- ((f1+f2+f3)*(NAD_3.NAM-NAM_3.NAM)+f4*(NAM_S3.NAM-NAM_3.NAM))/cNAM
    dNAM_4.NAM <- ((f1+f2+f3)*(0-NAM_4.NAM)+f4*(NAM_S4.NAM-NAM_4.NAM))/cNAM
    ####NA tracer####    
    dNAM_S0.NA <- -0.002
    dNAM_S6.NA <- 0.002
    dNA_0.NA <- 0
    dNA_6.NA <- 0
    dTrp_0.NA <- 0
    dTrp_6.NA <- 0
    dNAD_0.NA <- (f1*(Trp_0.NA-NAD_0.NA)+f2*(NA_0.NA-NAD_0.NA)+
                    f3*(NAM_0.NA-NAD_0.NA))/cNAD
    dNAD_6.NA <- (f1*(Trp_6.NA-NAD_6.NA)+
                    f2*(NA_6.NA-NAD_6.NA)+
                    f3*(NAM_6.NA-NAD_6.NA))/cNAD
    dNAM_0.NA <- ((f1+f2+f3)*(NAD_0.NA-NAM_0.NA)+f4*(NAM_S0.NA-NAM_0.NA))/cNAM
    dNAM_6.NA <- ((f1+f2+f3)*(NAD_6.NA-NAM_6.NA)+f4*(NAM_S6.NA-NAM_6.NA))/cNAM
    ####Tryptophan tracer#### 
    dNAM_S0.Trp <- 0.0018*Time^2-0.012*Time+0.0022
    dNAM_S6.Trp <- -0.0018*Time^2+0.012*Time-0.0022
    dTrp_0.Trp <- 0
    dTrp_6.Trp <- 0
    dNA_0.Trp <- 0
    dNA_6.Trp <- 0
    dNAD_0.Trp <- (f1*(Trp_0.Trp-NAD_0.Trp)+f2*(NA_0.Trp-NAD_0.Trp)+
                     f3*(NAM_0.Trp-NAD_0.Trp))/cNAD
    dNAD_6.Trp <- (f1*(Trp_6.Trp-NAD_6.Trp)+f2*(NA_6.Trp-NAD_6.Trp)+
                     f3*(NAM_6.Trp-NAD_6.Trp))/cNAD
    dNAM_0.Trp <- ((f1+f2+f3)*(NAD_0.Trp-NAM_0.Trp)+f4*(NAM_S0.Trp-NAM_0.Trp))/cNAM
    dNAM_6.Trp <- ((f1+f2+f3)*(NAD_6.Trp-NAM_6.Trp)+f4*(NAM_S6.Trp-NAM_6.Trp))/cNAM
    return(list(c(dNAM_S0.NAM,dNAM_S3.NAM,dNAM_S4.NAM,
                  dNA_0.NAM,dNA_3.NAM,dNA_4.NAM,
                  dTrp_0.NAM,dTrp_3.NAM,dTrp_4.NAM,
                  dNAD_0.NAM,dNAD_3.NAM,dNAM_0.NAM,dNAM_3.NAM,dNAM_4.NAM,
                  dNAM_S0.NA,dNAM_S6.NA,dNA_0.NA,dNA_6.NA,dTrp_0.NA,dTrp_6.NA,
                  dNAD_0.NA,dNAD_6.NA,dNAM_0.NA,dNAM_6.NA,
                  dNAM_S0.Trp,dNAM_S6.Trp,dNA_0.Trp,dNA_6.Trp,dTrp_0.Trp,dTrp_6.Trp,
                  dNAD_0.Trp,dNAD_6.Trp,dNAM_0.Trp,dNAM_6.Trp)))
  })
}


yini0 <- c(NAM_S0.NAM=1,NAM_S3.NAM=0,NAM_S4.NAM=0,NA_0.NAM=1,NA_3.NAM=0,NA_4.NAM=0,Trp_0.NAM=1,Trp_3.NAM=0,Trp_4.NAM=0,
           NAD_0.NAM=1,NAD_3.NAM=0,NAM_0.NAM=1,NAM_3.NAM=0,NAM_4.NAM=0,
           NAM_S0.NA=1,NAM_S6.NA=0,NA_0.NA=0.15,NA_6.NA=0.85,Trp_0.NA=1,Trp_6.NA=0,
           NAD_0.NA=1,NAD_6.NA=0,NAM_0.NA=1,NAM_6.NA=0,
           NAM_S0.Trp=1,NAM_S6.Trp=0,NA_0.Trp=1,NA_6.Trp=0,Trp_0.Trp=0.40,Trp_6.Trp=0.60,
           NAD_0.Trp=1,NAD_6.Trp=0,NAM_0.Trp=1,NAM_6.Trp=0)

times <- seq(0, 5, by = 0.01)

NADNAM <- function(f) {
  pars <- c(f1=f[1],f2=f[2],f3=f[3],f4=f[4],cNAD=NAD.Concentration,cNAM=NAM.Concentration)
  NADNAM1 <- tail(ode(func = NADNAMmod0D1, y = yini0, parms=pars, times = times[1:51],method="lsoda"),1)
  yini1 <- NADNAM1[2:35]
  names(yini1) <- colnames(NADNAM1)[2:35]
  NADNAM2 <- ode(func = NADNAMmod0D2, y = yini1, parms=pars, times = times[51:501],method="lsoda")
  return(list(t(NADNAM2[c(51,151,451),13:15]),NADNAM2[451,c(25,35)],t(NADNAM2[c(51,151,451),11:12]),NADNAM2[451,c(23,33)]))
}


Res <- function(f) {
  Results <- unlist(NADNAM(f))
  NAM.Diff <- Results[1:11]-unlist(NAM.Measurement)
  NAM.Diff <- NAM.Diff[c(2,3,5,6,8:11)]
  NAD.Diff <- Results[12:19]-unlist(NAD.Measurement)
  NAD.Diff <- NAD.Diff[c(2,3,5:8)]
  
  sum(NAM.Diff^2)/(0.015^2)+
    sum(NAD.Diff^2)/(0.008^2)
}


#Tissue List
#Adipose Brain Heart Kidney  Liver	Lung Pancreas Muscle Small Intestine Spleen

InputData <- read.csv("D:/XSu/Ling_NAD_Paper/Final Release/Figure5_Input.csv",header = FALSE,stringsAsFactors = FALSE)
FluxMatrix <- matrix(0,ncol=10,nrow=5)
colnames(FluxMatrix) <- c("Adipose","Brain","Heart","Kidney","Lung","Muscle","Pancreas","Small Intestine","Spleen","Liver")
rownames(FluxMatrix) <- c("f1","f2","f3","f4","res")

system.time(for (TissueCounter in 1:10) {
  
  NAM.Concentration <- as.numeric(as.character(InputData[TissueCounter*8-6,6]))
  NAD.Concentration <- as.numeric(as.character(InputData[TissueCounter*8-5,6]))
  
  NAM.Measurement <- list(ExpA=as.matrix(InputData[(TissueCounter*8-6):(TissueCounter*8-4),2:4]),
                          ExpB=as.numeric(InputData[TissueCounter*8-2,c(7,6)]))
  NAD.Measurement <- list(ExpA=as.matrix(InputData[(TissueCounter*8-3):(TissueCounter*8-2),2:4]),
                          ExpB=as.numeric(InputData[TissueCounter*8-1,c(7,6)]))
  
  system.time(fluxes<-bobyqa(c(10,10,300,50),lower=c(0,0,0,0),upper=c(100,20,800,100),fn=Res))
  
  FluxMatrix[1:4,TissueCounter] <- fluxes$par
  FluxMatrix[5,TissueCounter] <- fluxes$fval
  
  write.csv(FluxMatrix,file="Figure5_Output.csv")
})

fluxes$fval+qchisq(c(0.8,0.9,0.95),1)

FluxPath <- matrix(0,nrow=201,ncol=5)
FluxPath[1,1:4] <- fluxes$par
FluxPath[1,5] <- Res(fluxes$par)

CurrentFlux <- fluxes$par

system.time(for (i in 1:200) {
  J <- grad(Res,CurrentFlux)
  H <- hessian(Res,CurrentFlux)
  h <- 1
  b <- -J[c(1,2,4)]-H[c(1,2,4),3]*h
  A <- H[c(1,2,4),c(1,2,4)]
  u <- ginv(A,tol=1e-5) %*% b
  update <- c(u[1:2],h,u[3])
  CurrentFlux <- CurrentFlux + update
  FluxPath[1+i,1:4] <- CurrentFlux
  FluxPath[1+i,5] <- Res(CurrentFlux)
  flush.console()
  if(i%%1==0) print(paste(FluxPath[1+i,]))
  if(FluxPath[1+i,5]>20) break
})
write.csv(FluxPath,file="Save_Liver_High.csv")

FluxPath <- matrix(0,nrow=101,ncol=5)
FluxPath[1,1:4] <- fluxes$par
FluxPath[1,5] <- Res(fluxes$par)

CurrentFlux <- fluxes$par

system.time(for (i in 1:100) {
  J <- grad(Res,CurrentFlux)
  H <- hessian(Res,CurrentFlux)
  h <- -1
  b <- -J[c(1,2,4)]-H[c(1,2,4),3]*h
  A <- H[c(1,2,4),c(1,2,4)]
  u <- ginv(A,tol=1e-5) %*% b
  update <- c(u[1:2],h,u[3])
  CurrentFlux <- CurrentFlux + update
  FluxPath[1+i,1:4] <- CurrentFlux
  FluxPath[1+i,5] <- Res(CurrentFlux)
  if(FluxPath[1+i,5]>20) break
  flush.console()
  if(i%%1==0) print(paste(FluxPath[1+i,]))
})
write.csv(FluxPath,file="Save_Liver_Low.csv")