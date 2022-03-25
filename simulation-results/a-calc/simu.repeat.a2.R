##
## CALCULATE a AT FIXED P
##

{
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)
  
  library(mixmeta)
  library(survival)
  
  files.sources <- list.files(path = "funs/")
  x <- sapply(paste0("funs/", files.sources), source)
  N <- 30; tK0 <- 2; b0 <- 1; p <- 0.6
  S0 <- 100
}


ncores <- detectCores()
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2021)

##
## Scen 1; S1----
##
x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2


# S0 <- round(N/p)
DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0, b = b0, Unif.min = 50, Unif.max = 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_1 <- mean(DATA)

 
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0,  b = b0, Unif.min = 50, Unif.max = 150, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_2 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0,  b = b0, Unif.min = 100, Unif.max = 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_3 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0,  b = b0, Unif.min = 100, Unif.max = 300, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_4 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0, b = b0, Unif.min = 20, Unif.max = 60, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_5 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0,  b = b0, Unif.min = 20, Unif.max = 60, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_6 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0, b = b0, Unif.min = 300, Unif.max = 500, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_7 <- mean(DATA)


DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.p(S = S0, tK = tK0,  b = b0, Unif.min = 300, Unif.max = 500, CD = "EXP", Exp.r = 0.5, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p)
  
}
s1.a_8 <- mean(DATA)

parallel::stopCluster(cl)


a <- round(c(a1 = s1.a_1, a2 = s1.a_2, a3 = s1.a_3, a4 = s1.a_4, a5 = s1.a_5, a6 = s1.a_6, a7 = s1.a_7, a8 = s1.a_8),2)

save(a, file = paste0(b0, "_a.RData"))
