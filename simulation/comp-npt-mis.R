##******************************************************************************
##
## SIMULATION 
## WHEN CENSORING DISTRIBUTION IS INCORRECTLY SPECIFIED
## P=0.7 OR 0.5
##
##******************************************************************************

library("foreach")
library("parallel")
library("doSNOW")
library("doRNG")
library("survival")
library("magrittr")
library("numDeriv")

source("llk.TNM.ml.R")
source("llk.BNM.ml.R")
source("sroc.R")
source("sim.ipd.pdata.R")


simu.ipd2 <- function(S.n, npt.n=1, marker.n, censd=c("LN","UNIF")){ 
  
  data = .ps.sim.logit(
    S=set.ls$S[S.n],
    tK=set.ls$tK0,
    npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
    cens.dist = censd,  ## CENSOR DIST.
    marker.par=set.ls$marker[[marker.n]],
    alpha.n=marker.n,
    beta=5
  )
  
  pdata = data$pdata
  sdata = data$sdata
  pp = nrow(sdata)/nrow(pdata)
  
  fit3 <- BNM.mle(pdata$y1, pdata$y2, pdata$v1, pdata$v2, pdata$v12)
  mle3 <- c(fit3$par, beta=NA, conv = fit3$convergence, sauc = fit3$sauc*100, snp=nrow(pdata))
  
  fit4 <- BNM.mle(sdata$y1, sdata$y2, sdata$v1, sdata$v2, sdata$v12)
  mle4 <- c(fit4$par, beta=NA, conv = fit4$convergence, sauc = fit4$sauc*100, snp=nrow(sdata))
  
  initial.value = c(fit4$par,round(runif(1,0.5,6.5),1))
  initial.value[is.na(initial.value)] = round(runif(sum(is.na(initial.value)), 0.1,0.5),1)
  
  fit5 <- TNM.mle.sa(sdata$y1, sdata$y2, sdata$y3, sdata$v1, sdata$v2, sdata$v3, sdata$v12, sdata$v13, sdata$v23, pp, initial.values=initial.value)
  mle5 <- c(fit5$par, conv = fit5$convergence, sauc = fit5$sauc*100, snp=pp)
  

  ## COMBINE RES ------------------------------------------------
  
  x   <- cbind(mle3,mle4, mle5)
  rownames(x) <- c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3","beta","conv", "sauc", "snp")
  colnames(x) <- c("bnmp","bnmo","tnmsa")
  
  return(x)
  
}  

##
## 1000 REPEAT
##

ncores <- detectCores()-1
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

##p=0.7
for (npt.n in 1:2){
  load(paste0("nptmis-LN/alpha-npt",npt.n,".RData"))
  
  set.seed(2021)
  for(Sn in 1:3){  ## i list  # Sample Size
    for(marker.no in 1:5){
      
      DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("survival", "magrittr","numDeriv"), .errorhandling="remove")  %dorng%  {
        
        simu.ipd2(Sn, npt.n, marker.no, censd="LN")
        
      }
      save(DATA, file = paste0("nptmis-LN/simipd_npt", npt.n, "_marker",marker.no ,"_S",Sn,".RData"))
    }}}

##p=0.5
for (npt.n in 1:2){
  load(paste0("nptmis-UNIF/alpha-npt",npt.n,".RData"))
  
  set.seed(2021)
  for(Sn in 1:3){  ## i list  # Sample Size
    for(marker.no in 1:5){
      
      DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("survival", "magrittr","numDeriv"), .errorhandling="remove")  %dorng%  {
        
        simu.ipd2(Sn, npt.n, marker.no, censd="UNIF")
        
      }
      save(DATA, file = paste0("nptmis-UNIF/simipd_npt", npt.n, "_marker",marker.no ,"_S",Sn,".RData"))
    }}}

parallel::stopCluster(cl)






