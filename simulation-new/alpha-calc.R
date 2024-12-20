##
## COMPARE BNM SROC AND TNM SROC
##
rm(list = ls())
library("foreach")
library("parallel")
library("doSNOW")
library("doRNG")
library("survival")
library("magrittr")
library("numDeriv")

source("TNM.mle.R")
source("BNM.mle.R")
source("sroc.R")
source("sim.ipd.pdata.R")

beta=5

##******************************************************************************
##
## ALPHA WHEN CONSORING DISTRIBUTION IS CORRECTLY SPECIFIED
##
##******************************************************************************
ncores <- detectCores()-1
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)

## alpha for p=0.7
for(npt.n in 1:2){
  
  alpha = NULL
  for(marker.n in 1:6){
    
    
    set.seed(2021)
    a <- foreach(r=1:1000, .combine = "c", .packages = c("survival", "magrittr"))  %dorng%  {
      
      pdata <- .p.sim.logit(
        S=100,
        tK=set.ls$tK0,
        npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
        cens.dist = "EXP",  ## CENSOR DIST.
        marker.par = set.ls$marker[[marker.n]])
      
      pdata$t <- with(pdata, (y3/sqrt(v3)))
      
      fa = function(x) mean(pnorm(x + beta*pdata$t), na.rm=T)-0.7
      uniroot(fa, c(-10,5), extendInt="yes")$root
      
    }
    
    alpha = c(alpha, mean(a, na.rm = T))
    
  }
  
  save(alpha, file = paste0("7npt/alpha-npt",npt.n,".RData"))
}

## alpha for p=0.5
for(npt.n in 1:2){
  
  alpha = NULL
  for(marker.n in 1:6){
    
    
    set.seed(2021)
    a <- foreach(r=1:1000, .combine = "c", .packages = c("survival", "magrittr"))  %dorng%  {
      
      pdata <- .p.sim.logit(
        S=100,
        tK=set.ls$tK0,
        npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
        cens.dist = "EXP",  ## CENSOR DIST.
        marker.par = set.ls$marker[[marker.n]])
      
      pdata$t <- with(pdata, (y3/sqrt(v3)))
      
      fa = function(x) mean(pnorm(x + beta*pdata$t), na.rm=T)-0.5
      uniroot(fa, c(-10,5), extendInt="yes")$root
      
    }
    
    alpha = c(alpha, mean(a, na.rm = T))
    
  }
  
  save(alpha, file = paste0("5npt/alpha-npt",npt.n,".RData"))
}




##******************************************************************************
##
## ALPHA WHEN CONSORING DISTRIBUTION IS INCORRECTLY SPECIFIED
##
##******************************************************************************



## cens.dist = "LN"
for(npt.n in 1:2){
  
  alpha = NULL
  for(marker.n in 1:6){
    
    
    set.seed(2021)
    a <- foreach(r=1:1000, .combine = "c", .packages = c("survival", "magrittr"))  %dorng%  {
      
      pdata <- .p.sim.logit(
        S=100,
        tK=set.ls$tK0,
        npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
        cens.dist = "LN",  ## CENSOR DIST.
        marker.par = set.ls$marker[[marker.n]])
      
      pdata$t <- with(pdata, (y3/sqrt(v3)))
      
      fa = function(x) mean(pnorm(x + beta*pdata$t), na.rm=T)-0.7
      uniroot(fa, c(-10,5), extendInt="yes")$root
      
    }
    
    alpha = c(alpha, mean(a, na.rm = T))
    
  }
  
  save(alpha, file = paste0("nptmis-LN/alpha-npt",npt.n,".RData"))
}

## cens.dist = "UNIF"
for(npt.n in 1:2){
  
  alpha = NULL
  for(marker.n in 1:6){
    
    
    set.seed(2021)
    a <- foreach(r=1:1000, .combine = "c", .packages = c("survival", "magrittr"))  %dorng%  {
      
      pdata <- .p.sim.logit(
        S=100,
        tK=set.ls$tK0,
        npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
        cens.dist = "UNIF",  ## CENSOR DIST.
        marker.par = set.ls$marker[[marker.n]])
      
      pdata$t <- with(pdata, (y3/sqrt(v3)))
      
      fa = function(x) mean(pnorm(x + beta*pdata$t), na.rm=T)-0.7
      uniroot(fa, c(-10,5), extendInt="yes")$root
      
    }
    
    alpha = c(alpha, mean(a, na.rm = T))
    
  }
  
  save(alpha, file = paste0("nptmis-UNIF/alpha-npt",npt.n,".RData"))
}

parallel::stopCluster(cl)
