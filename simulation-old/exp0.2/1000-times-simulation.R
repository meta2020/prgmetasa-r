##******************************************************************************
##
## SIMULATION TO COMPARE COVERAGE % OF TNM IN POPULATION
##
##******************************************************************************

{
  setwd("~/GitHub-Bios/prgmetasa-r/simulation/exp0.2")
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)
  library(survival)
  library(testit)
  library(ucminf)
  
  # library(mixmeta)
  ## LOAD FUNCTIONS
  
  files.sources <- list.files(path = "../R/")
  x <- sapply(paste0("../R/", files.sources), source)
  rep <- 1

}

tK0 <- 2

dist1 <- c(N.min = 50, N.max = 150)
dist2 <- c(N.min = 50, N.max = 300)
pts.dist <- dist1

b0 <- 1



ncores <- 16*7
cl <- makeCluster(ncores, type="SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2023)


##******************************************************************************
##
## DISTRIBUTION OF PATIENTS: 1 ----
##
##******************************************************************************

pts.dist <- dist1
ptdis <- "ptdis1"
load(sprintf("a-calc/beta%s_pts%s_%s.RData", b0, pts.dist[1], pts.dist[2]))

for (S0 in c(200, 70)){
  
  ## BIOMARKER 1----
  alpha <- set$alpha1
  mk <- set$mk1
  hr <- "HR1"
  DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
    
    simu(S = S0, tK = tK0, b = b0, a = alpha, N.min = pts.dist[1], N.max= pts.dist[2], 
         CD = "EXP", Exp.r = 0.2, 
         x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
         true.p = NULL)
    
  }
  save(DATA, file = sprintf("S%d/%s_%s.RData", S0, ptdis, hr))
  
  ## BIOMARKER 2----
  alpha <- (set$alpha2)
  mk <- set$mk2
  hr <- "HR2"
  DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
    
    simu(S = S0, tK = tK0, b = b0, a = alpha, N.min = pts.dist[1], N.max= pts.dist[2], 
         CD = "EXP", Exp.r = 0.2, 
         x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
         true.p = NULL)
    
  }
  save(DATA, file = sprintf("S%d/%s_%s.RData", S0, ptdis, hr))


}


##******************************************************************************
##
## DISTRIBUTION OF PATIENTS: 2 ----
##
##******************************************************************************

pts.dist <- dist2
ptdis <- "ptdis2"
load(sprintf("a-calc/beta%s_pts%s_%s.RData", b0, pts.dist[1], pts.dist[2]))


for (S0 in c(200, 70)){
  
  ## BIOMARKER 1----
  alpha <- set$alpha1
  mk <- set$mk1
  hr <- "HR1"
  DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
    
    simu(S = S0, tK = tK0, b = b0, a = alpha, N.min = pts.dist[1], N.max= pts.dist[2], 
         CD = "EXP", Exp.r = 0.2, 
         x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
         true.p = NULL)
    
  }
  save(DATA, file = sprintf("S%d/%s_%s.RData", S0, ptdis, hr))
  
  ## BIOMARKER 2----
  alpha <- (set$alpha2)
  mk <- set$mk2
  hr <- "HR2"
  DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
    
    simu(S = S0, tK = tK0, b = b0, a = alpha, N.min = pts.dist[1], N.max= pts.dist[2], 
         CD = "EXP", Exp.r = 0.2, 
         x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
         true.p = NULL)
    
  }
  save(DATA, file = sprintf("S%d/%s_%s.RData", S0, ptdis, hr))
  
  
}


parallel::stopCluster(cl)






