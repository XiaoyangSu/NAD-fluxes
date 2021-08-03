require(deSolve)
require(minqa)
require(numDeriv)
require(MASS)

setwd("E:/XSu/Papers/NAD and Aging/2020May/Jul11")

NADNAMmod0D1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ####NAM tracer####
    dNAM_S0.NAM <- 0.752^Time*log(0.752)
    dNAM_S3.NAM <- 0
    dNAM_S4.NAM <- -0.752^Time*log(0.752)
    dTrp_0.NAM <- 0
    dTrp_3.NAM <- 0
    dTrp_4.NAM <- 0
    dNAD_0.NAM <- (f1*(Trp_0.NAM-NAD_0.NAM)+
                     f3*(NAM_0.NAM-NAD_0.NAM))/cNAD
    dNAD_3.NAM <- (f1*(Trp_3.NAM+Trp_4.NAM-NAD_3.NAM)+
                     f3*(NAM_3.NAM+NAM_4.NAM-NAD_3.NAM))/cNAD
    dNAM_0.NAM <- ((f1+f3)*(NAD_0.NAM-NAM_0.NAM)+f2*(NAM_S0.NAM-NAM_0.NAM))/cNAM
    dNAM_3.NAM <- ((f1+f3)*(NAD_3.NAM-NAM_3.NAM)+f2*(NAM_S3.NAM-NAM_3.NAM))/cNAM
    dNAM_4.NAM <- ((f1+f3)*(0-NAM_4.NAM)+f2*(NAM_S4.NAM-NAM_4.NAM))/cNAM
    
    ####Tryptophan tracer#### 
    dNAM_S0.Trp <- -0.021
    dNAM_S6.Trp <- 0.021
    dTrp_0.Trp <- 0
    dTrp_6.Trp <- 0
    dNAD_0.Trp <- (f1*(Trp_0.Trp-NAD_0.Trp)+
                     f3*(NAM_0.Trp-NAD_0.Trp))/cNAD
    dNAD_6.Trp <- (f1*(Trp_6.Trp-NAD_6.Trp)+
                     f3*(NAM_6.Trp-NAD_6.Trp))/cNAD
    dNAM_0.Trp <- ((f1+f3)*(NAD_0.Trp-NAM_0.Trp)+f2*(NAM_S0.Trp-NAM_0.Trp))/cNAM
    dNAM_6.Trp <- ((f1+f3)*(NAD_6.Trp-NAM_6.Trp)+f2*(NAM_S6.Trp-NAM_6.Trp))/cNAM
    return(list(c(dNAM_S0.NAM,dNAM_S3.NAM,dNAM_S4.NAM,
                  dTrp_0.NAM,dTrp_3.NAM,dTrp_4.NAM,
                  dNAD_0.NAM,dNAD_3.NAM,dNAM_0.NAM,dNAM_3.NAM,dNAM_4.NAM,
                  dNAM_S0.Trp,dNAM_S6.Trp,dTrp_0.Trp,dTrp_6.Trp,
                  dNAD_0.Trp,dNAD_6.Trp,dNAM_0.Trp,dNAM_6.Trp)))
  })
}

NADNAMmod0D2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ####NAM tracer####
    dNAM_S0.NAM <- -0.0774/Time
    dNAM_S3.NAM <- 0.0392/Time
    dNAM_S4.NAM <- 0.0382/Time
    dTrp_0.NAM <- 0
    dTrp_3.NAM <- 0
    dTrp_4.NAM <- 0
    dNAD_0.NAM <- (f1*(Trp_0.NAM-NAD_0.NAM)+
                     f3*(NAM_0.NAM-NAD_0.NAM))/cNAD
    dNAD_3.NAM <- (f1*(Trp_3.NAM+Trp_4.NAM-NAD_3.NAM)+
                     f3*(NAM_3.NAM+NAM_4.NAM-NAD_3.NAM))/cNAD
    dNAM_0.NAM <- ((f1+f3)*(NAD_0.NAM-NAM_0.NAM)+f2*(NAM_S0.NAM-NAM_0.NAM))/cNAM
    dNAM_3.NAM <- ((f1+f3)*(NAD_3.NAM-NAM_3.NAM)+f2*(NAM_S3.NAM-NAM_3.NAM))/cNAM
    dNAM_4.NAM <- ((f1+f3)*(0-NAM_4.NAM)+f2*(NAM_S4.NAM-NAM_4.NAM))/cNAM
    
    ####Tryptophan tracer#### 
    dNAM_S0.Trp <- -0.021
    dNAM_S6.Trp <- 0.021
    dTrp_0.Trp <- 0
    dTrp_6.Trp <- 0
    dNAD_0.Trp <- (f1*(Trp_0.Trp-NAD_0.Trp)+
                     f3*(NAM_0.Trp-NAD_0.Trp))/cNAD
    dNAD_6.Trp <- (f1*(Trp_6.Trp-NAD_6.Trp)+
                     f3*(NAM_6.Trp-NAD_6.Trp))/cNAD
    dNAM_0.Trp <- ((f1+f3)*(NAD_0.Trp-NAM_0.Trp)+f2*(NAM_S0.Trp-NAM_0.Trp))/cNAM
    dNAM_6.Trp <- ((f1+f3)*(NAD_6.Trp-NAM_6.Trp)+f2*(NAM_S6.Trp-NAM_6.Trp))/cNAM
    return(list(c(dNAM_S0.NAM,dNAM_S3.NAM,dNAM_S4.NAM,
                  dTrp_0.NAM,dTrp_3.NAM,dTrp_4.NAM,
                  dNAD_0.NAM,dNAD_3.NAM,dNAM_0.NAM,dNAM_3.NAM,dNAM_4.NAM,
                  dNAM_S0.Trp,dNAM_S6.Trp,dTrp_0.Trp,dTrp_6.Trp,
                  dNAD_0.Trp,dNAD_6.Trp,dNAM_0.Trp,dNAM_6.Trp)))
  })
}

NADNAMmod0D3 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ####NAM tracer####
    dNAM_S0.NAM <- -0.0774/Time
    dNAM_S3.NAM <- 0.0392/Time
    dNAM_S4.NAM <- 0.0382/Time
    dTrp_0.NAM <- 0
    dTrp_3.NAM <- 0
    dTrp_4.NAM <- 0
    dNAD_0.NAM <- (f1*(Trp_0.NAM-NAD_0.NAM)+
                     f3*(NAM_0.NAM-NAD_0.NAM))/cNAD
    dNAD_3.NAM <- (f1*(Trp_3.NAM+Trp_4.NAM-NAD_3.NAM)+
                     f3*(NAM_3.NAM+NAM_4.NAM-NAD_3.NAM))/cNAD
    dNAM_0.NAM <- ((f1+f3)*(NAD_0.NAM-NAM_0.NAM)+f2*(NAM_S0.NAM-NAM_0.NAM))/cNAM
    dNAM_3.NAM <- ((f1+f3)*(NAD_3.NAM-NAM_3.NAM)+f2*(NAM_S3.NAM-NAM_3.NAM))/cNAM
    dNAM_4.NAM <- ((f1+f3)*(0-NAM_4.NAM)+f2*(NAM_S4.NAM-NAM_4.NAM))/cNAM
    
    ####Tryptophan tracer#### 
    dNAM_S0.Trp <- -0.0478/Time
    dNAM_S6.Trp <- 0.0478/Time
    dTrp_0.Trp <- 0
    dTrp_6.Trp <- 0
    dNAD_0.Trp <- (f1*(Trp_0.Trp-NAD_0.Trp)+
                     f3*(NAM_0.Trp-NAD_0.Trp))/cNAD
    dNAD_6.Trp <- (f1*(Trp_6.Trp-NAD_6.Trp)+
                     f3*(NAM_6.Trp-NAD_6.Trp))/cNAD
    dNAM_0.Trp <- ((f1+f3)*(NAD_0.Trp-NAM_0.Trp)+f2*(NAM_S0.Trp-NAM_0.Trp))/cNAM
    dNAM_6.Trp <- ((f1+f3)*(NAD_6.Trp-NAM_6.Trp)+f2*(NAM_S6.Trp-NAM_6.Trp))/cNAM
    return(list(c(dNAM_S0.NAM,dNAM_S3.NAM,dNAM_S4.NAM,
                  dTrp_0.NAM,dTrp_3.NAM,dTrp_4.NAM,
                  dNAD_0.NAM,dNAD_3.NAM,dNAM_0.NAM,dNAM_3.NAM,dNAM_4.NAM,
                  dNAM_S0.Trp,dNAM_S6.Trp,dTrp_0.Trp,dTrp_6.Trp,
                  dNAD_0.Trp,dNAD_6.Trp,dNAM_0.Trp,dNAM_6.Trp)))
  })
}





yini0 <- c(NAM_S0.NAM=1,NAM_S3.NAM=0,NAM_S4.NAM=0,Trp_0.NAM=1,Trp_3.NAM=0,Trp_4.NAM=0,
           NAD_0.NAM=1,NAD_3.NAM=0,NAM_0.NAM=1,NAM_3.NAM=0,NAM_4.NAM=0,
           NAM_S0.Trp=1,NAM_S6.Trp=0,Trp_0.Trp=0.41,Trp_6.Trp=0.59,
           NAD_0.Trp=1,NAD_6.Trp=0,NAM_0.Trp=1,NAM_6.Trp=0)

times <- seq(0, 24, by = 0.01)

NADNAM <- function(f) {
  pars <- c(f1=f[1],f2=f[2],f3=f[3],cNAD=NAD.Concentration,cNAM=NAM.Concentration)
  NADNAM1 <- ode(func = NADNAMmod0D1, y = yini0, parms=pars, times = times[1:101],method="lsoda")[c(1,101),]
  yini1 <- NADNAM1[2,2:20]
  names(yini1) <- colnames(NADNAM1)[2:20]
  NADNAM2 <- ode(func = NADNAMmod0D2, y = yini1, parms=pars, times = times[101:401],method="lsoda")[seq(101,301,by=100),]
  yini2 <- NADNAM2[3,2:20]
  names(yini2) <- colnames(NADNAM2)[2:20]
  NADNAM3 <- ode(func = NADNAMmod0D3, y = yini2, parms=pars, times = times[401:2401],method="lsoda")[seq(101,2001,by=100),]
  return(rbind(NADNAM1,NADNAM2,NADNAM3))
  #return(list(t(NADNAM2[c(51,151,451),13:15]),NADNAM2[451,c(25,35)],t(NADNAM2[c(51,151,451),11:12]),NADNAM2[451,c(23,33)]))
}


Res <- function(f) {
  Results <- NADNAM(f)
  NAM.Diff <- Results[c(3,9,25),c(10:12,8,9)]-NAM.Measurement
  Trp.Diff <- Results[16,c(20,18)]-Trp.Measurement
  return(sum((NAM.Diff^2)/NAM.SD^2)+sum((Trp.Diff^2)/Trp.SD^2))
}

NADNAM2 <- function(f) {
  pars <- c(f1=f[1],f2=f[2],f3=f[3],cNAD=NAD.Concentration,cNAM=NAM.Concentration)
  NADNAM1 <- ode(func = NADNAMmod0D1, y = yini0, parms=pars, times = times[1:101],method="lsoda")
  yini1 <- NADNAM1[101,2:20]
  names(yini1) <- colnames(NADNAM1)[2:20]
  NADNAM2 <- ode(func = NADNAMmod0D2, y = yini1, parms=pars, times = times[101:401],method="lsoda")
  yini2 <- NADNAM2[301,2:20]
  names(yini2) <- colnames(NADNAM2)[2:20]
  NADNAM3 <- ode(func = NADNAMmod0D3, y = yini2, parms=pars, times = times[401:2401],method="lsoda")
  Result.Matrix <- rbind(NADNAM1,NADNAM2,NADNAM3)
  return(Result.Matrix)
  #return(list(t(NADNAM2[c(51,151,451),13:15]),NADNAM2[451,c(25,35)],t(NADNAM2[c(51,151,451),11:12]),NADNAM2[451,c(23,33)]))
}


InputData <- read.csv("E:/XSu/Papers/NAD and Aging/2020May/Jul11/Input_SDAvg_Young_0711.csv",header = FALSE,stringsAsFactors = FALSE)
FluxMatrix <- matrix(0,ncol=12,nrow=4)
colnames(FluxMatrix) <- c("Liver","Kidney","Quad","BAT","Heart","Spleen","Jejunum","Brain","Pancreas","iWAT","Ileum","PC")
rownames(FluxMatrix) <- c("f1","f2","f3","res")

NAM.Measurement <- matrix(0,ncol=5,nrow=3)
NAM.SD <- matrix(0,ncol=5,nrow=3)


system.time(for (TissueCounter in 1:12) {
  
  NAM.Concentration <- as.numeric((InputData[TissueCounter*7-1,7]))
  NAD.Concentration <- as.numeric((InputData[TissueCounter*7,7]))
  
  NAM.Measurement <- apply(as.matrix(InputData[(TissueCounter*7-c(5,3,1)),2:6]), 2, as.numeric)
  NAM.SD <- apply(as.matrix(InputData[(TissueCounter*7-c(4,2,0)),2:6]), 2, as.numeric)
  
  
  Trp.Measurement <-as.numeric((InputData[TissueCounter*7-c(5,3),7]))
  Trp.SD <-as.numeric((InputData[TissueCounter*7-c(4,2),7]))
  
  fluxes <- DEoptim(fn = Res,lower=rep(0,3),upper=c(100,200,5000),control = list(itermax=200, parallelType=1,
                parVar=c("NADNAM","NADNAMmod0D1","NADNAMmod0D2","NADNAMmod0D3","NAD.Concentration","NAM.Concentration",
                               "times","NAM.Measurement","NAM.SD","Trp.Measurement","Trp.SD","yini0","ode")))
  
  FluxMatrix[1:3,TissueCounter] <- fluxes$optim$bestmem
  FluxMatrix[4,TissueCounter] <- fluxes$optim$bestval
  
  Fluxes.Results <- fluxes$optim$bestmem
  names(Fluxes.Results) <- NULL
  
  write.csv(NADNAM2(Fluxes.Results),file=paste(colnames(FluxMatrix)[TissueCounter],"_Young_SDAvg_TimeCourse.csv",sep=""))
  write.csv(FluxMatrix,file="Output_Young_SDAvg_2020Jul11.csv")
})


##########The actual code ENDS here. The rest is doodle code##################


##############################################################################
system.time(fluxes<-bobyqa(c(10,10,300),lower=c(0,0,0),upper=c(100,200,800),fn=Res))



write.csv(NADNAM2(fluxes$par),file="Kinetics.csv")

Flux.NoSD <- fluxes


########Test code
pars <- c(f1=0.001,f2=0.002,f3=0.003,cNAD=NAD.Concentration,cNAM=NAM.Concentration)
NADNAM1 <- ode(func = NADNAMmod0D1, y = yini0, parms=pars, times = times[1:101],method="lsoda")

N <- NADNAM(c(0.001,0.002,0.003))


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