
for (p in seq(0.9, 0.1, -0.1)){
S0 <- round(N/p)
DATA <- foreach(r = 1:1000, .combine = "cbind", .packages = c("survival"), .errorhandling = "remove")  %dorng%  {  # 
  
  simu(S = S0, tK = tK0, b = b0, a = alpha[p*10], N.min = pts.dist[1], N.max= pts.dist[2], 
       CD = "UNIF", minc=1, maxc=4, 
       x1.u = mk[1], x2.u = mk[2], x1.s = mk[3], x2.s = mk[4], v.sd = mk[5],
       true.p = NULL)
  
}
save(DATA, file = sprintf("N%d/%s_%s_p%.1f.RData", N, ptdis, hr, p))

}
