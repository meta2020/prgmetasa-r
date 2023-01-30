
# 0.9 ----
p0 <- 0.9
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0, b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.9 <- mean(DATA)

# 0.8 ----
p0 <- 0.8
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0, b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.8 <- mean(DATA)

# 0.7 ----
p0 <- 0.7
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "c", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0, b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.7 <- mean(DATA)

# 0.6 ----
p0 <- 0.6 
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.6 <- mean(DATA)

# 0.5 ----
p0 <- 0.5 
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.5 <- mean(DATA)

# 0.4 ----
p0 <- 0.4 
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0,N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.4 <- mean(DATA)

# 0.3 ----
p0 <- 0.3 
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0,N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.3 <- mean(DATA)

# 0.2 ----
p0 <- 0.2
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.2 <- mean(DATA)

# 0.1 ----
p0 <- 0.1
S0 <- round(N/p0)
DATA <- foreach(r=1:1000, .combine = "cbind", .packages = c("mixmeta", "survival"))  %dorng%  {  # 
  
  sim.alpha(S = S0, tK = tK0,  b = b0, N.min = pts.dist[1], N.max = pts.dist[2], CD = "UNIF", minc=1, maxc=4, x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5], p = p0)
  
}
hr0.1 <- mean(DATA)