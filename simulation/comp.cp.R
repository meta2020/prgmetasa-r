##
## COMPARE BNM SROC AND TNM SROC
##
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


##
## SIMULATION SCENARIOS
##

set.ls <- list(
  S = c(35, 50, 200),
  tK0 = 2,
  npt = list(c(50,100),
             c(50,200)),
  marker = list(c(x1.u = 0.7, x2.u = 0.3, x1.s = 0.1, x2.s = 0.3, v.sd = 0.14), ## HZ
                c(x1.u = 0.7, x2.u = 0.3, x1.s = 0.1, x2.s = 0.3, v.sd = 0.2), ## change var
                c(x1.u = 0.7, x2.u = 0.3, x1.s = 0.7, x2.s = 0.9, v.sd = 0.2), ##
                c(x1.u = 1, x2.u = 0.2, x1.s = 0.7, x2.s = 0.9, v.sd = 0.2), ## SMALL HR
                c(x1.u = 1, x2.u = 0.2, x1.s = 0.7, x2.s = 0.9, v.sd = 0.3)) ## SMALL HR
)


##
## CALCULATE ALPHA 
##
# 
# ncores <- detectCores()-1
# cl <- makeCluster(ncores, "SOCK")
# doSNOW::registerDoSNOW(cl)
# 
# 
# alpha = NULL
# for(marker.n in 1:5){
# 
#   set.seed(2021)
#   a <- foreach(r=1:1000, .combine = "c", .packages = c("survival", "magrittr"))  %dorng%  {
# 
#     pdata <- .p.sim.logit(
#         S=100,
#         tK=set.ls$tK0,
#         npt.range=set.ls$npt[[2]],  ## NUMBER OF PTS
#         cens.dist = "EXP",  ## CENSOR DIST.
#         marker.par = set.ls$marker[[marker.n]])
# 
#   pdata$t <- with(pdata, y3/sqrt(v3))
# 
#   fa = function(x) mean(pnorm(x + 3*pdata$t), na.rm=T)-0.7
#   uniroot(fa, c(-10,5), extendInt="yes")$root
# 
#   }
# 
#   alpha = c(alpha, mean(a, na.rm = T))
# 
# }
# 
# parallel::stopCluster(cl)
# #
# save(alpha, file = "alpha-npt2.RData")


##
## SIMULATION FUNCTION
##
simu.ipd <- function(S.n, npt.n=1, marker.n, beta=3){ 
  
  data = .ps.sim.logit(
    S=set.ls$S[S.n],
    tK=set.ls$tK0,
    npt.range=set.ls$npt[[npt.n]],  ## NUMBER OF PTS
    cens.dist = "EXP",  ## CENSOR DIST.
    marker.par=set.ls$marker[[marker.n]],
    alpha.n=marker.n,
    beta
  )
  
  pdata = data$pdata
  sdata = data$sdata
  pp = nrow(sdata)/nrow(pdata)
  
  fit3 <- BNM.mle(pdata$y1, pdata$y2, pdata$v1, pdata$v2, pdata$v12)
  mle3 <- c(fit3$par, beta=NA, conv = fit3$convergence, sauc = fit3$sauc*100, sauc.lb=NA, sauc.ub=NA, snp=nrow(pdata), cp=NA)
  
  fit4 <- BNM.mle(sdata$y1, sdata$y2, sdata$v1, sdata$v2, sdata$v12)
  mle4 <- c(fit4$par, beta=NA, conv = fit4$convergence, sauc = fit4$sauc*100, sauc.lb=fit4$sauc.lb*100, sauc.ub=fit4$sauc.ub*100, snp=nrow(sdata), 
            cp = fit3$sauc*100 >=fit4$sauc.lb*100 && fit3$sauc*100 <=fit4$sauc.ub*100)
  
  initial.value = c(fit4$par,NA)
  initial.value[is.na(initial.value)] = runif(sum(is.na(initial.value)), -0.1,0.1)
  
  fit5 <- TNM.mle.sa(sdata$y1, sdata$y2, sdata$y3, sdata$v1, sdata$v2, sdata$v3, sdata$v12, sdata$v13, sdata$v23, pp, initial.values=initial.value)
  mle5 <- c(fit5$par, conv = fit5$convergence, sauc = fit5$sauc*100, sauc.lb=fit5$sauc.lb*100, sauc.ub=fit5$sauc.ub*100, snp=pp, 
            cp = fit3$sauc*100 >=fit5$sauc.lb*100 && fit3$sauc*100 <=fit5$sauc.ub*100)
  

  ## COMBINE RES ------------------------------------------------
  
  x   <- cbind(mle3,mle4, mle5)
  rownames(x) <- c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3","beta","conv", "sauc", "sauc.lb","sauc.ub","snp","cp")
  colnames(x) <- c("bnmp","bnmo","tnmsa")
  
  return(x)
  
}  

##
## 1000 REPEAT
##
load("alpha-npt1.RData")
ncores <- detectCores()-1
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)


set.seed(2021)

for(Sn in 1:3){  ## i list  # Sample Size
  for(marker.no in 1){
    ## 1
    DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("survival", "magrittr","numDeriv"))  %dorng%  {
      
      simu.ipd(Sn, npt.n=1, marker.no, beta=3)
      # simu.ipd(1,1,1,3)
      
    }
    save(DATA, file = paste0("simipd_npt1_marker",marker.no ,"_S",Sn,".RData"))
  }}

parallel::stopCluster(cl)


##
## SUMMARY
##
Sn=1
mark.n=1

.table.pars2 <- function(Sn, mark.n){
  
  load(paste0("simipd_npt1_marker",mark.n ,"_S",Sn,".RData"))
  DATA[c(1:10, 12:14,16), DATA[11,]!=0] <- NA
  dim1 <- nrow(DATA)
  dim2 <- ncol(DATA)/1000
  
  name.c <- colnames(DATA)[1:dim2]
  
  dim(DATA) <- c(dim1, dim2, 1000)
  
  med <- apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))
  sd  <- apply(DATA, 1:2, function(x) sd(x, na.rm = TRUE))
  rownames(med) <- c("mu1","mu2","mu3", "tau1","tau2","tau3", "rho1","rho2","rho3", "beta","conv", "sauc","sauc.lb","sauc.ub","SNp","cp")
  med[11,] = 1-med[11,]
  # med[12,] = med[12,]-med[12,1]
  
  return(round(med[c(12,15,11,16),],2))
}

cbind(
  .table.pars2(1,1),
  .table.pars2(2,1),
  .table.pars2(3,1))
cbind(
  .table.pars2(1,2),
  .table.pars2(2,2),
  .table.pars2(3,2))
cbind(
  .table.pars2(1,3),
  .table.pars2(2,3),
  .table.pars2(3,3))
cbind(
  .table.pars2(1,4),
  .table.pars2(2,4),
  .table.pars2(3,4))
cbind(
  .table.pars2(1,5),
  .table.pars2(2,5),
  .table.pars2(3,5))

