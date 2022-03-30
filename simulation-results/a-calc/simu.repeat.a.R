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
  N <- 100; tK0 <- 2; b0 <- 1
}


ncores <- detectCores()
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2021)

##******************************************************************************
##
## Scen 1; S1---- MODERATE HR
##
##******************************************************************************

x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0, b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s1.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s1.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0,N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s1.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s1.a_0.2 <- mean(DATA)


##******************************************************************************
##
## Scen 1; S2---- SMALL HR
##
##******************************************************************************

x1.u <- 0.6; x2.u <- 0.4; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s2.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s2.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4 
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s2.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)
  
}
s2.a_0.2 <- mean(DATA)



##******************************************************************************
##
## Scen 1; S3---- LARGE HR
##
##******************************************************************************

x1.u <- 0.8; x2.u <- 0.2; x1.s <- 0.3; x2.s <- 0.1; v.sd <- 0.2

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
s3.a_0.8 <- mean(DATA)

# 0.6 ----
p0 <- 0.6
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
s3.a_0.6 <- mean(DATA)

# 0.4 ----
p0 <- 0.4
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
s3.a_0.4 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:100, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  #

  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = 50, N.max= 150, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd, p = p0)

}
s3.a_0.2 <- mean(DATA)



parallel::stopCluster(cl)



a <- round(
  c(s1.a_0.8,s1.a_0.6,s1.a_0.4,s1.a_0.2,
    s2.a_0.8,s2.a_0.6,s2.a_0.4,s2.a_0.2,
    s3.a_0.8,s3.a_0.6,s3.a_0.4,s3.a_0.2),2)

names(a) <- c(paste0("hr1.", seq(0.8, 0.2,-0.2)), paste0("hr2.", seq(0.8, 0.2,-0.2)), paste0("hr3.", seq(0.8, 0.2,-0.2)))

save(a, file = paste0("beta_", b0, "_alpha_", N, ".RData"))
