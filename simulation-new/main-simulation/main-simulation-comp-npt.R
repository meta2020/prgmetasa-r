##******************************************************************************
##
## SIMULATION 
## WHEN CENSORING DISTRIBUTION IS CORRECTLY SPECIFIED
## P=0.7 OR 0.5
##
##******************************************************************************
rm(list = ls())

simu.ipd <- function(S.n, npt.n, marker.n, p = p){ 
  
  data = .ps.sim.logit(
    S=set.ls$S[S.n],
    tK=set.ls$tK0,
    npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
    cens.dist = "EXP",  ## CENSOR DIST.
    marker.par=set.ls$marker[[marker.n]],
    alpha.n=marker.n,
    beta=5
  )
  
  pdata = data$pdata
  sdata = data$sdata
  pp = nrow(sdata)/nrow(pdata)
  
  fit3 <- BNM.mle(pdata$y1, pdata$y2, pdata$v1, pdata$v2, pdata$v12)
  mle3 <- c(fit3$par, beta=NA, conv = fit3$convergence, sauc = fit3$sauc, snp=nrow(pdata))
  
  fit1 <- TNM.mle(pdata$y1, pdata$y2, pdata$y3, pdata$v1, pdata$v2, pdata$v3, pdata$v12, pdata$v13, pdata$v23)
  mle1 <- c(fit1$par, beta=NA, conv = fit1$convergence, sauc = fit1$sauc, snp=nrow(pdata))
  
  fit4 <- BNM.mle(sdata$y1, sdata$y2, sdata$v1, sdata$v2, sdata$v12)
  mle4 <- c(fit4$par, beta=NA, conv = fit4$convergence, sauc = fit4$sauc, snp=nrow(sdata))
  
  fit2 <- TNM.mle(sdata$y1, sdata$y2, sdata$y3, sdata$v1, sdata$v2, sdata$v3, sdata$v12, sdata$v13, sdata$v23)
  mle2 <- c(fit2$par, beta=NA, conv = fit2$convergence, sauc = fit2$sauc, snp=nrow(sdata))
  
  initial.value = c(fit4$par,round(runif(1,4.5,6.5),1))
  initial.value[is.na(initial.value)] = round(runif(sum(is.na(initial.value)), 0.1,0.5),1)
  
  fit5 <- TNM.mle.sa(sdata$y1, sdata$y2, sdata$y3, sdata$v1, sdata$v2, sdata$v3, sdata$v12, sdata$v13, sdata$v23, pp, initial.values=initial.value)
  mle5 <- c(fit5$par, conv = fit5$convergence, sauc = fit5$sauc, snp=pp)
  

  ## COMBINE RES ------------------------------------------------
  
  x   <- cbind(mle3,mle1,mle4,mle2, mle5)
  rownames(x) <- c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3","beta","conv", "sauc", "snp")
  colnames(x) <- c("bnmp","tnmp","bnmo","tnmo","tnmsa")
  
  return(x)
  
}  

library("foreach")
library("parallel")
library("doSNOW")
library("doRNG")
library("survival")
library("magrittr")
library("numDeriv")

source("BNM.mle.R")
source("TNM.mle.R")
source("TNM.mle.sa.R")
source("sroc.R")
source("sim.ipd.pdata.R")

##
## 1000 REPEAT
##

ncores <- detectCores()-1
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

##p=0.7
for (npt.n in 1:2){
load(paste0("../7npt/alpha-npt",npt.n,".RData"))
  
set.seed(2021)
for(Sn in 1:3){  ## i list  # Sample Size
  for(marker.no in 1:6){

    DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("survival", "magrittr","numDeriv"), .errorhandling="remove")  %dorng%  {
      
      simu.ipd(Sn, npt.n, marker.no, p=0.7)

    }
    save(DATA, file = paste0("../7npt/simipd_npt", npt.n, "_marker",marker.no ,"_S",Sn,".RData"))
  }}}

##p=0.5
for (npt.n in 1:2){
  load(paste0("../5npt/alpha-npt",npt.n,".RData"))
  
  set.seed(2021)
  for(Sn in 4:6){  ## i list  # Sample Size
    for(marker.no in 1:6){
      
      DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("survival", "magrittr","numDeriv"), .errorhandling="remove")  %dorng%  {
        simu.ipd(Sn, npt.n, marker.no, p=0.5)
        
      }
      save(DATA, file = paste0("../5npt/simipd_npt", npt.n, "_marker",marker.no ,"_S",Sn,".RData"))
    }}
  }

parallel::stopCluster(cl)



