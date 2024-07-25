##******************************************************************************
##
## SIMULATION 
## WHEN CENSORING DISTRIBUTION IS CORRECTLY SPECIFIED
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

set.seed(123)
S.n=3
npt.n=1
marker.n=5

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
mle3 <- c(fit3$par, beta=NA, conv = fit3$convergence, sauc = fit3$sauc*100, snp=nrow(pdata))

fit4 <- BNM.mle(sdata$y1, sdata$y2, sdata$v1, sdata$v2, sdata$v12)
mle4 <- c(fit4$par, beta=NA, conv = fit4$convergence, sauc = fit4$sauc*100, snp=nrow(sdata))

initial.value = c(fit4$par,round(runif(1,0.5,6.5),1))
initial.value[is.na(initial.value)] = round(runif(sum(is.na(initial.value)), 0.1,0.5),1)

fit5 <- TNM.mle.sa(sdata$y1, sdata$y2, sdata$y3, sdata$v1, sdata$v2, sdata$v3, sdata$v12, sdata$v13, sdata$v23, pp, initial.values=initial.value)
mle5 <- c(fit5$par, conv = fit5$convergence, sauc = fit5$sauc*100, snp=pp)
x   <- cbind(mle3,mle4, mle5)
rownames(x) <- c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3","beta","conv", "sauc", "snp")
colnames(x) <- c("bnmp","bnmo","tnmsa")

plot(plogis(-pdata$y2), plogis(pdata$y1), xlim = c(0,1), ylim = c(0,1))
points(plogis(-sdata$y2), plogis(sdata$y1), xlim = c(0,1), ylim = c(0,1), col=2)
.SROC(par=x[c(1,2,4,5,7),], add = T, ncols=1:3, add.spoint=F)





