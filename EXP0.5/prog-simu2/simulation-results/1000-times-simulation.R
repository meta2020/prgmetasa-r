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
  load(sprintf("a-calc/beta_1_alpha_%d.RData", N))

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

x1.u <- set$set1[1]; x2.u <- set$set1[2]; x1.s <- set$set1[3]; x2.s <- set$set1[4]; v.sd <- set$set1[5]


# 0.8 ----
p <- 0.8
S0 <- round(N/p)
alpha <- set$alpha1[1]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  # 
	
	simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR1_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6 
S0 <- round(N/p)
alpha <- set$alpha1[2]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR1_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4 
S0 <- round(N/p)
alpha <- set$alpha1[3]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR1_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
alpha <- set$alpha1[4]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR1_p%.1f.RData", N, p))





##******************************************************************************
##
## Scen 1; S2----
##
##******************************************************************************


x1.u <- set$set2[1]; x2.u <- set$set2[2]; x1.s <- set$set2[3]; x2.s <- set$set2[4]; v.sd <- set$set2[5]

# 0.8 ----
p <- 0.8
S0 <- round(N/p)
alpha <- set$alpha2[1]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("N%d/HR2_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6
S0 <- round(N/p)
alpha <- set$alpha2[2]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("N%d/HR2_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4
S0 <- round(N/p)
alpha <- set$alpha2[3]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("N%d/HR2_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
alpha <- set$alpha2[4]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #

  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

}
save(DATA, file = sprintf("N%d/HR2_p%.1f.RData", N, p))





##******************************************************************************
##
## Scen 1; S3 ---- LARGE HR
##
##******************************************************************************

x1.u <- set$set3[1]; x2.u <- set$set3[2]; x1.s <- set$set3[3]; x2.s <- set$set3[4]; v.sd <- set$set3[5]

# 0.8 ----
p <- 0.8
S0 <- round(N/p)
alpha <- set$alpha3[1]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR3_p%.1f.RData", N, p))

# 0.6 ----
p <- 0.6
S0 <- round(N/p)
alpha <- set$alpha3[2]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR3_p%.1f.RData", N, p))

# 0.4 ----
p <- 0.4
S0 <- round(N/p)
alpha <- set$alpha3[3]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR3_p%.1f.RData", N, p))

# 0.2 ----
p <- 0.2
S0 <- round(N/p)
alpha <- set$alpha3[4]
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "remove")  %dorng%  {  #
  
  simu(S = S0, tK = tK0, b = 1, a = alpha, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("N%d/HR3_p%.1f.RData", N, p))


parallel::stopCluster(cl)





