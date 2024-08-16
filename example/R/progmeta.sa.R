##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE) USING EXPANSION
##
##******************************************************************************


llk.TNM.ml = function(par, y1, y2, y3, v1, v2, v3, v12, v13, v23){
  
  mu = par[1:3]  
  tau = par[4:6]     
  rho = par[7:9]    
  
  t11 = tau[1]^2
  t22 = tau[2]^2
  t33 = tau[3]^2
  
  t12 = prod(tau[-3])*rho[1]
  t23 = prod(tau[-1])*rho[2]
  t13 = prod(tau[-2])*rho[3]
  
  y=cbind(y1,y2,y3)
  v=cbind(v1,v2,v3,v12,v13,v23)
  
  f=rep(NA, nrow(y))
  for (i in 1:nrow(y)){
    mu.sigma=matrix(c(t11+v[i,1], t12+v[i,4], t13+v[i,5], 
                      t12+v[i,4], t22+v[i,2], t23+v[i,6],
                      t13+v[i,5], t23+v[i,6], t33+v[i,3]),3,3)
    if(is.null(pd.solve(mu.sigma, silent = T))) f[i]= NA else f[i] = dmnorm(y[i,], mu, mu.sigma, log=TRUE)  
  }
  
  s.l1=sum(f, na.rm = TRUE)
  
  return(-s.l1)
}


clk.TNM.ml = function(
    par, y1, y2, y3,
    v1, v2, v3, v12, v13, v23,
    p, a.interval
    
){
  
  mu = par[1:3]  
  tau = par[4:6]     
  rho = par[7:9]    
  
  
  t11 = tau[1]^2
  t22 = tau[2]^2
  t33 = tau[3]^2
  
  t12 = prod(tau[-3])*rho[1]
  t23 = prod(tau[-1])*rho[2]
  t13 = prod(tau[-2])*rho[3]
  
  beta = par[10] 
  
  t_lnHR   = abs(y3/sqrt(v3))
  
  ##
  ## FUNCTOIN b(Sigma) ----
  ##
  
  f.b = function(a){
    
    sq = suppressWarnings(sqrt(1 + beta^2*(1 + t33/v3)))
    pnorm( (a + beta * mu[3]/sqrt(v3)) / sq)
    
  }
  
  
  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##
  
  a.p = function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}
  
  a.opt.try = suppressWarnings(try(uniroot(a.p, interval=a.interval, extendInt="yes"), silent = TRUE)) 
  
  a.opt = a.opt.try$root
  
  
  f.l1=NULL
  y=cbind(y1,y2,y3)
  v=cbind(v1,v2,v3,v12,v13,v23)
  for (i in 1:nrow(y)){

    mu.sigma=matrix(c(t11+v[i,1], t12+v[i,4], t13+v[i,5], 
                      t12+v[i,4], t22+v[i,2], t23+v[i,6],
                      t13+v[i,5], t23+v[i,6], t33+v[i,3]),3,3)
    if(is.null(pd.solve(mu.sigma, silent = T))) f.l1[i]= NA else f.l1[i] = dmnorm(y[i,], mu, mu.sigma, log=TRUE)  
    }
  
  s.l1=sum(f.l1, na.rm = TRUE)
  
  ##
  ##  LOGLIKELIHOOD-2 OF a(a.opt) ----
  ##
  
  f.l2 = pnorm(a.opt + beta * t_lnHR, log.p = T)
  s.l2 = sum(f.l2, na.rm = TRUE)
  
  
  ##
  ##  LOGLIKELIHOOD-3 OF b(a.opt) ----
  ##
  
  f.l3 = f.b(a.opt)
  s.l3 = sum(log(f.l3), na.rm = TRUE)
  
  
  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##
  
  return(-(s.l1 + s.l2 - s.l3)) ## NEGATIVE
  
  
}

progmeta.sa = function(
    y1, y2, y3, v1, v2, v3, v12, v13, v23, 
    p, 
    a.interval= c(-5, 5)
    ){

  
  ## initial value

  fn.int = function(par) llk.TNM.ml(par, y1, y2, y3, v1, v2, v3, v12, v13, v23)
  opt.int = nlminb(c(rep(0.1,6), rep(-0.1, 3)), fn.int,
                  lower = c(rep(-Inf,2), -Inf, rep(0,3), rep(-1,3)),
                  upper = c(rep( Inf,2),  Inf, rep(Inf,3), rep( 1,3)))
  init.value = c(opt.int$par,0.5)
  
  ## optimization 

  fn.o = function(par) clk.TNM.ml(par, y1, y2, y3, v1, v2, v3, v12, v13, v23,
                                  p, a.interval)

  opt.o = nlminb(init.value, 
                 fn.o,
                 lower = c(rep(-Inf,3),rep(0,3), rep(-0.999,3), 0.001),
                 upper = c(rep( Inf,3),rep(Inf,3), rep(0.999,3), 5))

  
  num.hessian = hessian(fn.o, opt.o$par)
  rownames(num.hessian) = colnames(num.hessian) = c(paste0("mu",1:3), paste0("tau",1:3), paste0("rho",1:3), "beta")

  if(p==1) var.ml = solve(num.hessian[-10,-10]) else var.ml = solve(num.hessian)
  var.matrix = var.ml[c(1,2,4,5,7), c(1,2,4,5,7)]

  sauc.ci = SAUC.ci(opt.o$par[c(1,2,4,5,7)], var.matrix = var.matrix, sauc.type = "sroc")

  beta   = opt.o$par[10]
  t33 = (opt.o$par[6])^2
  mu3  = opt.o$par[3]

  f.b = function(a){
    
    sq = suppressWarnings(sqrt(1 + beta^2 * (1 + t33 / v3)))
    pnorm( (a + beta * mu3/sqrt(v3)) / sq )
    
  }
  
  a.p = function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}
  
  a.opt.try = suppressWarnings(try(uniroot(a.p, interval=a.interval, extendInt="yes"), silent = TRUE)) 
  if (!inherits(a.opt.try, "try-error")) a.opt = a.opt.try$root else a.opt = NA
  p.hat = 1/mean(1/f.b(a.opt), na.rm = TRUE)
  
  par = c(opt.o$par,  a.opt, p.hat, sauc.ci, plogis(opt.o$par[1:2]), opt.o$convergence)
  names(par) = c(paste0("mu",1:3), paste0("tau",1:3), paste0("rho",1:3), "beta", "alpha",
                 "p.hat","sauc", "sauc.lb", "sauc.ub","sen", "spe", "conv")
  
return(par)

}


