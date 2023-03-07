##******************************************************************************
##
## SIMULATION TO COMPARE COVERAGE % OF TNM IN POPULATION
##
##******************************************************************************

{
  
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)
  library(survival)
  library(testit)
  ## LOAD FUNCTIONS
  
  files.sources <- list.files(path = "../R/")
  x <- sapply(paste0("../R/", files.sources), source)
  rep <- 1000
  
}

tK0 <- 2

dist1 <- c(N.min = 50, N.max = 150)
dist2 <- c(N.min = 50, N.max = 300)
pts.dist <- dist1

b0 <- 1



ncores <- 120
cl <- makeCluster(ncores)
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

for (N in c(20,30,50,100)){
  
  ## BIOMARKER 1----
  alpha <- rev(set$alpha1)
  mk <- set$mk1
  hr <- "HR1"
  source("loop-hr-incorrect.R")
  
  ## BIOMARKER 2----
  alpha <- rev(set$alpha2)
  mk <- set$mk2
  hr <- "HR2"
  source("loop-hr-incorrect.R")

  
}


##******************************************************************************
##
## DISTRIBUTION OF PATIENTS: 2 ----
##
##******************************************************************************

pts.dist <- dist2
ptdis <- "ptdis2"
load(sprintf("a-calc/beta%s_pts%s_%s.RData", b0, pts.dist[1], pts.dist[2]))


for (N in c(20,30,50,100)){
  
  ## BIOMARKER 1----
  alpha <- rev(set$alpha1)
  mk <- set$mk1
  hr <- "HR1"
  source("loop-hr-incorrect.R")
  
  ## BIOMARKER 2----
  alpha <- rev(set$alpha2)
  mk <- set$mk2
  hr <- "HR2"
  source("loop-hr-incorrect.R")
  
  
}


parallel::stopCluster(cl)






