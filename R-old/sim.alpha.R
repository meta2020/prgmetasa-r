##
## CALCULATE a, GIVEN p = p0
##

# S <- 30; tK <- 3; b <- 0.5; a <- -0.7; N.min <- 50; N.max <- 100; CD = "EXP"; Exp.r <- 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; p <- 0.7

sim.alpha <- function(
  S, tK, b, N.min, N.max, 
  CD, Exp.r, meanlog, sdlog, minc, maxc, 
  x1.u, x2.u, x1.s, x2.s, v.sd, 
  p
  ){
	
  ## POPULATION DATA
  
  p_dataSt.all <- population.data.list(S, tK, N.min, N.max, CD, Exp.r, meanlog, sdlog, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd)
  p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))
  
  ## REMOVE NA
  
  p_dataSt <- p_dataSt[complete.cases(p_dataSt), ]
  
  ## ANALYTICAL POPULATION DATA 

  eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
  dataSt.p <- convert.dt(data = p_dataSt, tK = tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)
  
  ## REMOVE NA
  
  dataSt.p <- dataSt.p[complete.cases(dataSt.p), ]
  
  ## 2. CALCULATE a ----
  
  fa <- function(a) mean(pnorm(a + b*dataSt.p$t_lnHR)) - p
  uniroot(fa, c(-10, 10), extendInt="yes")$root

}
