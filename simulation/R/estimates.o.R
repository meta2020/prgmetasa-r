estimates.o <- function(
    dataSt.o = dataSt.o.7, p.hat =p.hat.7, N = N.7, pp=pp.7, pe=pe.7,
    # bnm.intoo= bnm.into, tnm.intoo=tnm.into, 
    tnmsa.intpp = tnmsa.intp,
    bnm.oo = bnm.o, tnm.oo = tnm.o, 
    # tnm.intt = tnm.int,
    dataSt.pp = dataSt.p
    # tnm.into=tnm.int
    ){
  
  eps <- 0
  ## 3.1 BNM.O.7 ----
  
  {
    fn.bnm.ml  <- function(par) llk.BNM.ml(
      par,
      dataSt.o$u_sen, dataSt.o$u_spe,
      dataSt.o$v_sen, dataSt.o$v_spe, dataSt.o$v_senspe)

    fit.bnm_o  <- nlminb(
      # bnm.intoo,
      c(rep(0.1,4),-0.1),
      fn.bnm.ml,
      lower = c(rep(-5,2), rep(0, 2),-1),
      upper = c(rep( 5,2), rep(3, 2), 1)
    )

    if(!inherits(fit.bnm_o, "try-error")) {

      conv  <- fit.bnm_o$convergence
      sauc  <- SAUC(par = fit.bnm_o$par)  # u1 u2 t1 t2 r1
      bnm.oo[c(1,2,4,5,7,10, 11, 14, 15, 16)]  <- c(fit.bnm_o$par, sauc, conv, N, pp, pe)
    }
    
    # fit.bnm_o <- try(mixmeta::mixmeta(cbind(dataSt.o$u_sen, dataSt.o$u_spe),S=cbind(dataSt.o$v_sen, dataSt.o$v_senspe, dataSt.o$v_spe), data=dataSt.o, method = "ml"), silent = TRUE)
    # 
    # if((!inherits(fit.bnm_o, "try-error")) ) {
    #   
    #   conv <- fit.bnm_o$converged  # 0 is success
    #   mu <- fit.bnm_o$coefficients
    #   v  <- fit.bnm_o$Psi
    #   tau_rho <- c(sqrt(v[c(1,4)]), v[3]/prod(sqrt(v[c(1,4)])))
    #   sauc <- SAUC(par = c(mu, tau_rho))  # u1 u2 t1 t2 r1
    #   bnm.oo[c(1,2,4,5,7,10,11)]  <- c(mu, tau_rho, sauc, conv-1)
    #   
    # } 
  }
  
  
  ## 3.2 TNM.O  ----
  
  {
    fn.tnm.ml  <- function(par) llk.TNM.ml(
      par, 
      dataSt.o$u_sen, dataSt.o$u_spe, dataSt.o$u_lnHR,
      dataSt.o$v_sen, dataSt.o$v_spe, dataSt.o$v_lnHR,
      dataSt.o$v_senspe, dataSt.o$v_senlnHR, dataSt.o$v_spelnHR)
    
    # initial values
    if(!inherits(fit.bnm_o, "try-error")) {
      
      if (fit.bnm_o$convergence==0) {
        
        tnm.into <- bnm.oo[1:9]
        tnm.into[c(3,6)] <- c(0.1,0.1)#runif(2, 0.1, 1)
        tnm.into[c(8,9)] <- c(-0.1,-0.1)#runif(2, -1, 1)
        
      } else tnm.into <- c(rep(0.1,6), rep(-0.1,3))} else tnm.into <- c(rep(0.1,6), rep(-0.1,3))
    
    fit.tnm_o  <- nlminb(
      tnm.into,
      fn.tnm.ml,
      lower = c(rep(-5,3), rep(eps, 3),rep(-1,3)),
      upper = c(rep( 5,3), rep(3, 3),rep(1,3))
    )
    
    if(!inherits(fit.tnm_o, "try-error")) {
      
      conv  <- fit.tnm_o$convergence
      sauc  <- SAUC(par = fit.tnm_o$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
      tnm.oo[1:11] <- c(fit.tnm_o$par, sauc, conv)
      
      if (conv==0) tnmsa.into <- c(fit.tnm_o$par, 0.5) else tnmsa.into <- c(tnm.into, 0.5)
      
    } else tnmsa.into <- c(tnm.into, 0.5)
    
  }
  
  ## 3 TNM.SA.7  ----
  {
    fn.tnm.ml  <- function(par) clk.TNM.ml(
      par, 
      dataSt.o$u_sen, dataSt.o$u_spe, dataSt.o$u_lnHR,
      dataSt.o$v_sen, dataSt.o$v_spe, dataSt.o$v_lnHR,
      dataSt.o$v_senspe, dataSt.o$v_senlnHR, dataSt.o$v_spelnHR, 
      p=p.hat, a.interval = c(-5,5))
    
    #tnm.intp with initial values: tnmsa.intp
    fit.tnm_sa_p <- nlminb(
      tnmsa.intpp, fn.tnm.ml,
      lower = c(rep(-5,3), rep(eps,3), rep(-1,3), eps),
      upper = c(rep( 5,3), rep(3, 3), rep( 1,3), 2) )
    
    if(!inherits(fit.tnm_sa_p, "try-error")) {
      
      conv  <- fit.tnm_sa_p$convergence
      sauc  <- SAUC(par = fit.tnm_sa_p$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
      
      ## CALCULATE alpha.p
      # b <- 1
      b   <- fit.tnm_sa_p$par[10]
      u3  <- fit.tnm_sa_p$par[3]
      t33 <- fit.tnm_sa_p$par[6]^2
      
      a.p <- function(a) {mean(1/f.b(a, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE) - 1/p.hat}
      a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-5,5), extendInt="downX"), silent = FALSE))
      if (!inherits(a.opt.try, "try-error")) alpha.p <- a.opt.try$root else alpha.p <- NA
      
      p.est.a  <- mean(pnorm(alpha.p+b*(dataSt.pp$t_lnHR)), na.rm = TRUE)
      p.est.b  <- mean(f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.pp$v_lnHR),   na.rm = TRUE)
      p.est    <- 1/mean(1/f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE)
      
      tnm.sa_p <- c(fit.tnm_sa_p$par[1:9], sauc, conv, b, alpha.p, NA, p.est, p.est.a, p.est.b)
      
    }  
  }
  
  {
    #tnm.into with inital values: tnmsa.into
    fit.tnm_sa_o <- nlminb(
      tnmsa.into, fn.tnm.ml,
      lower = c(rep(-5,3), rep(eps, 3), rep(-1,3), eps),
      upper = c(rep( 5,3), rep(3, 3), rep( 1,3), 2) )
    
    if(!inherits(fit.tnm_sa_o, "try-error")) {
      
      conv  <- fit.tnm_sa_o$convergence
      sauc  <- SAUC(par = fit.tnm_sa_o$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
      
      ## CALCULATE alpha.p
      # b <- 1
      b   <- fit.tnm_sa_o$par[10]
      u3  <- fit.tnm_sa_o$par[3]
      t33 <- fit.tnm_sa_o$par[6]^2
      
      a.p <- function(a) {mean(1/f.b(a, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE) - 1/p.hat}
      a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-5,5), extendInt="downX"), silent = FALSE))
      if (!inherits(a.opt.try, "try-error")) alpha.p <- a.opt.try$root else alpha.p <- NA
      
      p.est.a  <- mean(pnorm(alpha.p+b*(dataSt.pp$t_lnHR)), na.rm = TRUE)
      p.est.b  <- mean(f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.pp$v_lnHR),   na.rm = TRUE)
      p.est    <- 1/mean(1/f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE)
      
      tnm.sa_o <- c(fit.tnm_sa_o$par[1:9], sauc, conv, b, alpha.p, NA, p.est, p.est.a, p.est.b)
      
    }  
  }
  
  {  
    #tnm.int with initial values: vague initial values
    fit.tnm_sa_n <- nlminb(
      c(c(rep(0.1,6), rep(-0.1, 3)), 0.5), fn.tnm.ml,
      lower = c(rep(-5,3), rep(eps, 3), rep(-1,3), eps),
      upper = c(rep( 5,3), rep(3, 3), rep( 1,3), 2) )
    
    if(!inherits(fit.tnm_sa_n, "try-error")) {
      
      conv  <- fit.tnm_sa_n$convergence
      sauc  <- SAUC(par = fit.tnm_sa_n$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
      
      ## CALCULATE alpha.p
      # b <- 1
      b   <- fit.tnm_sa_n$par[10]
      u3  <- fit.tnm_sa_n$par[3]
      t33 <- fit.tnm_sa_n$par[6]^2
      
      a.p <- function(a) {mean(1/f.b(a, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE) - 1/p.hat}
      a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-5,5), extendInt="downX"), silent = FALSE))
      if (!inherits(a.opt.try, "try-error")) alpha.p <- a.opt.try$root else alpha.p <- NA
      
      p.est.a  <- mean(pnorm(alpha.p+b*(dataSt.pp$t_lnHR)), na.rm = TRUE)
      p.est.b  <- mean(f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.pp$v_lnHR),   na.rm = TRUE)
      p.est    <- 1/mean(1/f.b(alpha.p, b=b, u3=u3, t33=t33, v3=dataSt.o$v_lnHR), na.rm = TRUE)
      
      tnm.sa_n <- c(fit.tnm_sa_n$par[1:9], sauc, conv, b, alpha.p, NA, p.est, p.est.a, p.est.b)
      
    }  
  }
  
  res <- cbind(bnm.oo, tnm.oo, tnm.sa_o, tnm.sa_p, tnm.sa_n)
  colnames(res) <- c("BNM.O", "TNM.O", "SA.o", "SA.p", "SA.n")
  rownames(res) <- c(
    "mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3", "sauc", "conv",
    "beta", "alpha.p", "SN", "p.est", "p.est.a", "p.est.b")
  return(res)
  
}


