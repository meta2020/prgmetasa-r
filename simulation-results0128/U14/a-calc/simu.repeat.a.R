##******************************************************************************
##
## CALCULATE ALPHA AT DIFFERENT MARGINAL SELECTION PROBBAILITY (P)
##
##******************************************************************************

## LOAD FUNCTIONS

{
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)

  library(survival)
  library(testit)
  files.sources <- list.files(path = "../../../R/")
  x <- sapply(paste0("../../../R/", files.sources), source)
}

## SET SCENARIOS ----

tK0 <- 2
b0 <- 1

## DIST OF PATIENTS ~ U(N.min, N.max)

dist1 <- c(N.min = 50, N.max = 150)
dist2 <- c(N.min = 50, N.max = 300)

## BIOMARKES
 
mk1 <- c(x1.u = 0.7, x2.u = 0.3, x1.s = 0.1, x2.s = 0.3, v.sd = 0.1) ## SMALL HR
mk2 <- c(x1.u = 0.6, x2.u = 0.4, x1.s = 0.2, x2.s = 0.4, v.sd = 0.1) ## MODERATE HR
mk3 <- c(x1.u = 0.8, x2.u = 0.2, x1.s = 0.1, x2.s = 0.3, v.sd = 0.05) ## LARGE HR

## STUDY NUMBER
N <- 50

## START REPEAT

ncores <- 120
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2023)

##******************************************************************************
##
## N = 50 ----
##
##******************************************************************************

pts.dist <- dist1
mk <- mk1
source("repeat.p.unif.R")
alpha1 <- c(hr1.a_0.9 = hr0.9, hr1.a_0.8 = hr0.8, hr1.a_0.7 = hr0.7, hr1.a_0.6 = hr0.6,
            hr1.a_0.5 = hr0.5, hr1.a_0.4 = hr0.4, hr1.a_0.3 = hr0.3, hr1.a_0.2 = hr0.2,
            hr1.a_0.1 = hr0.1)

mk <- mk2
source("repeat.p.unif.R")
alpha2 <- c(hr2.a_0.9 = hr0.9, hr2.a_0.8 = hr0.8, hr2.a_0.7 = hr0.7, hr2.a_0.6 = hr0.6,
            hr2.a_0.5 = hr0.5, hr2.a_0.4 = hr0.4, hr2.a_0.3 = hr0.3, hr2.a_0.2 = hr0.2,
            hr2.a_0.1 = hr0.1)

mk <- mk3
source("repeat.p.unif.R")
alpha3 <- c(hr3.a_0.9 = hr0.9, hr3.a_0.8 = hr0.8, hr3.a_0.7 = hr0.7, hr3.a_0.6 = hr0.6,
            hr3.a_0.5 = hr0.5, hr3.a_0.4 = hr0.4, hr3.a_0.3 = hr0.3, hr3.a_0.2 = hr0.2,
            hr3.a_0.1 = hr0.1)

set <- list(
  N = N,
  pts.dist = pts.dist,
  mk1 = mk1, 
  mk2 = mk2, 
  mk3 = mk3,
  beta = b0, 
  alpha1 = round(alpha1,2),
  alpha2 = round(alpha2,2),
  alpha3 = round(alpha3,3)
)

save(set, file = sprintf("beta%s_pts%s_%s.RData", b0, pts.dist[1], pts.dist[2]))

##******************************************************************************
##
## N = 50 ----
##
##******************************************************************************

pts.dist <- dist2
mk <- mk1
source("repeat.p.unif.R")
alpha1 <- c(hr1.a_0.9 = hr0.9, hr1.a_0.8 = hr0.8, hr1.a_0.7 = hr0.7, hr1.a_0.6 = hr0.6,
            hr1.a_0.5 = hr0.5, hr1.a_0.4 = hr0.4, hr1.a_0.3 = hr0.3, hr1.a_0.2 = hr0.2,
            hr1.a_0.1 = hr0.1)

mk <- mk2
source("repeat.p.unif.R")
alpha2 <- c(hr2.a_0.9 = hr0.9, hr2.a_0.8 = hr0.8, hr2.a_0.7 = hr0.7, hr2.a_0.6 = hr0.6,
            hr2.a_0.5 = hr0.5, hr2.a_0.4 = hr0.4, hr2.a_0.3 = hr0.3, hr2.a_0.2 = hr0.2,
            hr2.a_0.1 = hr0.1)

mk <- mk3
source("repeat.p.unif.R")
alpha3 <- c(hr3.a_0.9 = hr0.9, hr3.a_0.8 = hr0.8, hr3.a_0.7 = hr0.7, hr3.a_0.6 = hr0.6,
            hr3.a_0.5 = hr0.5, hr3.a_0.4 = hr0.4, hr3.a_0.3 = hr0.3, hr3.a_0.2 = hr0.2,
            hr3.a_0.1 = hr0.1)

set <- list(
  N = N,
  pts.dist = pts.dist,
  mk1 = mk1, 
  mk2 = mk2, 
  mk3 = mk3,
  beta = b0, 
  alpha1 = round(alpha1,2),
  alpha2 = round(alpha2,2),
  alpha3 = round(alpha3,3)
)

save(set, file = sprintf("beta%s_pts%s_%s.RData", b0, pts.dist[1], pts.dist[2]))


parallel::stopCluster(cl)
