

for(p in seq(0.9, 0.1, -0.1)){
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = b0, a = alpha[9], N.min = pts.dist[1], N.max= pts.dist[2], 
       CD = "EXP", Exp.r = 0.2, 
       x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
       true.p = NULL)
  
}
save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
}

# # 0.8 ----
# p <- 0.8
# S0 <- round(N/p)
# DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[8], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# # 0.7 ----
# p <- 0.7 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[7], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# 
# # 0.6 ----
# p <- 0.6 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[6], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# # 0.5 ----
# p <- 0.5 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[5], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# 
# # 0.4 ----
# p <- 0.4 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[4], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# # 0.3 ----
# p <- 0.3 
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[3], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# 
# # 0.2 ----
# p <- 0.2
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[2], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
# # 0.1 ----
# p <- 0.1
# S0 <- round(N/p)
# DATA <- foreach(r = 1:rep, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
#   
#   simu(S = S0, tK = tK0, b = b0, a = alpha[1], N.min = pts.dist[1], N.max = pts.dist[2], 
#        CD = "EXP", Exp.r = 0.2, 
#        x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
#        true.p = NULL)
#   
# }
# save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))
# 
