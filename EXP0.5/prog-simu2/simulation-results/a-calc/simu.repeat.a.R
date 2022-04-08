##
## CALCULATE ALPHA AT DIFFERENT MARGINAL SELECTION PROBBAILITY (P)
##

{
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)
  
  library(mixmeta)
  library(survival)
  library(testit)
  files.sources <- list.files(path = "../../funs/")
  x <- sapply(paste0("../../funs/", files.sources), source)
  N <- 30; tK0 <- 2; b0 <- 1
}


ncores <- detectCores()
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2021)

##******************************************************************************
##
## Scen 1; S1---- SMALL HR
##
##******************************************************************************

x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.3; x2.s <- 0.3; v.sd <- 0.2
set1 <- c(x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0, b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr1.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr1.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0,N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr1.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr1.a_0.2 <- mean(DATA)



##******************************************************************************
##
## Scen 1; S2---- SMALL HR
##
##******************************************************************************

x1.u <- 0.6; x2.u <- 0.4; x1.s <- 0.5; x2.s <- 0.5; v.sd <- 0.2
set2 <- c(x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr2.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr2.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr2.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
hr2.a_0.2 <- mean(DATA)



##******************************************************************************
##
## Scen 1; S3---- SMALL HR
##
##******************************************************************************

x1.u <- 0.5; x2.u <- 0.4; x1.s <- 0.3; x2.s <- 0.3; v.sd <- 0.2
set3 <- c(x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
hr3.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
hr3.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
hr3.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
hr3.a_0.2 <- mean(DATA)


parallel::stopCluster(cl)



alpha1 <- c(hr1.a_0.8 = hr1.a_0.8, hr1.a_0.6 = hr1.a_0.6, hr1.a_0.4 = hr1.a_0.4, hr1.a_0.2 = hr1.a_0.2)
alpha2 <- c(hr2.a_0.8 = hr2.a_0.8, hr2.a_0.6 = hr2.a_0.6, hr2.a_0.4 = hr2.a_0.4, hr2.a_0.2 = hr2.a_0.2)
alpha3 <- c(hr3.a_0.8 = hr3.a_0.8, hr3.a_0.6 = hr3.a_0.6, hr3.a_0.4 = hr3.a_0.4, hr3.a_0.2 = hr3.a_0.2)


set <- list(set1 = set1, alpha1 = alpha1,
            set2 = set2, alpha2 = alpha2,
            set3 = set3, alpha3 = alpha3)


save(set, file = paste0("beta_", b0, "_alpha_", N, ".RData"))
