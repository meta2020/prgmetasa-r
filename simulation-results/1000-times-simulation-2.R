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
  library(mixmeta)
  library(survival)
  library(testit)
  ## LOAD FUNCTIONS
  
  files.sources <- list.files(path = "../funs/")
  x <- sapply(paste0("../funs/", files.sources), source)
  N <- 30; tK0 <- 2
  load(sprintf("a-calc/beta_1_alpha2_%d.RData", N))

}


ncores <- detectCores() 
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2021)

##******************************************************************************
##
## Scen 1; S1----
##
##******************************************************************************

x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.5


# 0.8 ----
p <- 0.8
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
	
	simu(S = S0, tK = tK0, b = 1, a = a[1], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR1_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6 
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = a[2], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR1_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4 
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = a[3], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR1_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = a[4], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR1_p%.1f.RData", N, p))


# parallel::stopCluster(cl)



##******************************************************************************
##
## Scen 1; S2----
##
##******************************************************************************


x1.u <- 0.6; x2.u <- 0.4; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.5

# 0.8 ----
p <- 0.8
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = a[5], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("2N%d/HR2_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = a[6], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("2N%d/HR2_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = a[7], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("2N%d/HR2_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = a[8], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("2N%d/HR2_p%.1f.RData", N, p))


# parallel::stopCluster(cl)



##******************************************************************************
##
## Scen 1; S3 ---- LARGE HR
##
##******************************************************************************

x1.u <- 0.8; x2.u <- 0.2; x1.s <- 0.3; x2.s <- 0.1; v.sd <- 0.5

# 0.8 ----
p <- 0.8
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = a[9], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR3_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = a[10], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR3_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = a[11], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR3_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = a[12], N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("2N%d/HR3_p%.1f.RData", N, p))


parallel::stopCluster(cl)




