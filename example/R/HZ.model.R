##******************************************************************************
##
## THE HZ MODEL 
##
##******************************************************************************


llk.BNM.ml = function(par,y1, y2, v1, v2, v12){

  mu = par[1:2]
  tau = par[3:4]
  rho = par[5]

  t11 = tau[1]^2
  t22 = tau[2]^2
  t12 = prod(tau)*rho

  f=NULL
  y=cbind(y1,y2)
  v=cbind(v1,v2,v12)
  for (i in 1:nrow(y)){
    mu.sigma=matrix(c(t11+v[i,1], t12+v[i,3],
                      t12+v[i,3], t22+v[i,2]),2,2)
    if(is.null(pd.solve(mu.sigma, silent = T))) f[i]= NA else f[i] = dmnorm(y[i,], mu, mu.sigma, log=TRUE)  
  }
  
  s.l1=sum(f, na.rm = TRUE)
  return(-s.l1)
}
  
HZ.model = function(y1,y2,v1,v2,v12){

  fn = function(par) llk.BNM.ml(par, y1=y1, y2=y2, v1=v1, v2=v2, v12=v12)
  opt.bnm.ml = nlminb(c(rep(0.1,4), -0.1), fn, 
                       lower = c(rep(-Inf,2), rep(0,2), -0.999),
                       upper = c(rep( Inf,2), rep(Inf,2), 0.999))
  
  hes.bnm    <- numDeriv::hessian(fn, opt.bnm.ml$par)
  rownames(hes.bnm) <- colnames(hes.bnm) <- c("mu.se", "mu.sp", "tau.sp", "tau.sp", "rho")
  sauc.bnm   <- SAUC.ci(opt.bnm.ml$par, var.matrix = solve(hes.bnm), sauc.type = "sroc")
  
  par = c(opt.bnm.ml$par, sauc.bnm, plogis(opt.bnm.ml$par[1:2]), opt.bnm.ml$convergence)
  names(par) = c(paste0("mu",1:2), paste0("tau",1:2), "rho", "sauc", "sauc.lb", "sauc.ub","sen", "spe", "conv")
  
  return(par)
}
