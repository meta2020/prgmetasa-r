##------------------------------------------------------------
##
## SIMULATION TO COMPARE COVERAGE % OF TNM IN POPULATION
##
##------------------------------------------------------------

# setwd("~/Documents/meeting2/20211018/PB-simulation")

{
  
  rm(list = ls())
  
  library(foreach)
  library(parallel)
  library(doSNOW)
  library(doRNG)
  library(mixmeta)
  library(survival)
  
  files.sources <- list.files(path = "../funs/")
  x <- sapply(paste0("../funs/", files.sources), source)
  N <- 30; tK0 <- 2

}

load("../a-calc/1.5_a2_30.RData")
ncores <- detectCores() 
cl <- makeCluster(ncores, "SOCK")
doSNOW::registerDoSNOW(cl)
set.seed(2021)

##
## Scen 1; S1----
##

x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2


# 0.9 ----
p <- 0.9
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
	
	simu.run.PB(S = S0, tK = tK0, b = 1, a = a[1], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("resb2/DATA_N%d_HR1_p%.1f.RData", N, p))

# 0.7 ----
p <- 0.7 
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu.run.PB(S = S0, tK = tK0, b = 1, a = a[2], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("resb2/DATA_N%d_HR1_p%.1f.RData", N, p))

# 0.5 ----
p <- 0.5 
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu.run.PB(S = S0, tK = tK0, b = 1, a = a[3], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("resb2/DATA_N%d_HR1_p%.1f.RData", N, p))

# 0.3 ----
p <- 0.3
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu.run.PB(S = S0, tK = tK0, b = 1, a = a[4], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("resb2/DATA_N%d_HR1_p%.1f.RData", N, p))

# 0.1 ----
p <- 0.1
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
  
  simu.run.PB(S = S0, tK = tK0, b = 1, a = a[5], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
  
}
save(DATA, file = sprintf("resb2/DATA_N%d_HR1_p%.1f.RData", N, p))

parallel::stopCluster(cl)


# 
# 
# ##
# ## Scen 1; S2----
# ##
# 
# x1.u <- 0.6; x2.u <- 0.4; x1.s <- 0.3; x2.s <- 0.3; v.sd <- 0.2
# 
# # 0.9 ----
# p <- 0.9
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
#   
#   simu.run.PB(S = S0, tK = tK0, b = 1, a = a[6], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
#   
# }
# save(DATA, file = sprintf("resb2/DATA_N%d_HR2_p%.1f.RData", N, p))
# 
# # 0.7 ----
# p <- 0.7 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
#   
#   simu.run.PB(S = S0, tK = tK0, b = 1, a = a[7], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
#   
# }
# save(DATA, file = sprintf("resb2/DATA_N%d_HR2_p%.1f.RData", N, p))
# 
# # 0.5 ----
# p <- 0.5 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
#   
#   simu.run.PB(S = S0, tK = tK0, b = 1, a = a[8], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
#   
# }
# save(DATA, file = sprintf("resb2/DATA_N%d_HR2_p%.1f.RData", N, p))
# 
# # 0.3 ----
# p <- 0.3
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
#   
#   simu.run.PB(S = S0, tK = tK0, b = 1, a = a[9], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
#   
# }
# save(DATA, file = sprintf("resb2/DATA_N%d_HR2_p%.1f.RData", N, p))
# 
# # 0.1 ----
# p <- 0.1
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"), .errorhandling = "stop")  %dorng%  {  # 
#   
#   simu.run.PB(S = S0, tK = tK0, b = 1, a = a[10], Unif.min = 100, Unif.max= 300, CD = "EXP", Exp.r = 0.2, x1.u = x1.u, x2.u = x2.u, x1.s = x1.s, x2.s = x2.s, v.sd = v.sd)
#   
# }
# save(DATA, file = sprintf("resb2/DATA_N%d_HR2_p%.1f.RData", N, p))
# 
# 
# 
# 
# 
# 
