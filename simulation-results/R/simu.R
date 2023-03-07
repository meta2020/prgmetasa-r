##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu <- function(
    S, tK, b, a, N.min, N.max, CD, Exp.r, 
    meanlog, sdlog, minc, maxc, 
    x1.u, x2.u, x1.s, x2.s, v.sd, 
    true.p=NULL){

  
  # S = 50; tK = tK0; b = b0; a = set$alpha1[3]; N.min = pts.dist[1]; N.max= pts.dist[2]
  # CD = "EXP"; Exp.r = 0.2
  # x1.u = set$mk1[1]; x2.u = set$mk1[2]; x1.s = set$mk1[3]; x2.s = set$mk1[4]; v.sd = set$mk1[5]
  # true.p = NULL
  # eps <- sqrt(.Machine$double.eps)
  
   eps <- 0
  
  ####### RUN
  {
    ## POPULATION DATA
    
    p_dataSt.all <- population.data.list(
      S, tK, N.min, N.max, 
      CD, Exp.r, meanlog, sdlog, minc, maxc, 
      x1.u, x2.u, x1.s, x2.s, v.sd)
    p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))
    
    ## REMOVE NA
    
    p_dataSt <- na.omit(p_dataSt)
    
    ## ANALYTICAL POPULATION DATA 
    
    eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
    dataSt.p <- convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)
    
    ## REMOVE NA
    
    dataSt.p <- na.omit(dataSt.p)
    
    ## 2. CREATE OBSERVED LIST ----
    
    # b <- 3; a <- -1
    select.p  <- pnorm(b * (dataSt.p$t_lnHR) + a)
    dataSt.p$select.p  <- select.p
    dataSt.p$select.id <- sapply(1:nrow(dataSt.p), function(i){rbinom(1, 1, select.p[i]) })
    
    dataSt.o <- dataSt.p[dataSt.p$select.id == 1, ]
    
    
    SS <- nrow(dataSt.p)
    NN <- nrow(dataSt.o)
    
    ## MARGINAL p
    # true.p <- NULL
    
    prop <- nrow(dataSt.o)/nrow(dataSt.p) 
    p.e <- mean(pnorm(b * dataSt.p$t_lnHR + a), na.rm = TRUE)
    if (is.null(true.p)) p.hat <- p.e else p.hat <- true.p 
    
    bnm.p <- c(rep(NA, 11), SS, NN, prop, NA, NA)
    bnm.o <- c(rep(NA, 11), NA, NA, p.hat, NA, NA)
    tnm.sa_p <- tnm.sa_o <-tnm.sa <- rep(NA, 16)

    
    if(nrow(dataSt.o)==0) {
      
      res <- NULL
      
    } else {
    
    
    {
      ## POPULATION DATA
      
      p.y1  <- dataSt.p$u_sen
      p.y2  <- dataSt.p$u_spe 
      p.y3  <- dataSt.p$u_lnHR
      
      p.s11  <- dataSt.p$v_sen
      p.s22  <- dataSt.p$v_spe 
      p.s33  <- dataSt.p$v_lnHR
      
      p.s12 <- dataSt.p$v_senspe 
      p.s13 <- dataSt.p$v_senlnHR 
      p.s23 <- dataSt.p$v_spelnHR 
      t.p   <- dataSt.p$t_lnHR

      ## OBSERVED DATA
      
      o.y1  <- dataSt.o$u_sen
      o.y2  <- dataSt.o$u_spe 
      o.y3  <- dataSt.o$u_lnHR
      
      o.s11  <- dataSt.o$v_sen
      o.s22  <- dataSt.o$v_spe 
      o.s33  <- dataSt.o$v_lnHR
      
      o.s12 <- dataSt.o$v_senspe 
      o.s13 <- dataSt.o$v_senlnHR 
      o.s23 <- dataSt.o$v_spelnHR 
      
    }

    
    ## SET INTITIAL VALUES
    # bnm.int <- c(mean(o.y1), mean(o.y2),
    #              rep(0.1,2), -0.1)
    # tnm.int <- c(mean(o.y1), mean(o.y2), mean(o.y3),
    #              rep(0.1,3), rep(-0.1, 3))
                
    tnm.int <- c(rep(0,3), rep(0.1, 3), rep(-0.1, 3))
    bnm.int <- c(rep(0,2), rep(0.1, 2), -0.1)
    b0 <- 0.5
    
    
    ## 1 BNM.P ----
    {
      
      fn.bnm.ml  <- function(par) llk.BNM.ml(par, p.y1, p.y2, p.s11, p.s22, p.s12)
      fit.bnm_p  <- nlminb(bnm.int, fn.bnm.ml,
                           lower = c(rep(-Inf,2), rep(0, 2),  -1),
                           upper = c(rep( Inf,2), rep(Inf, 2), 0) 
                           )

      if(!inherits(fit.bnm_p, "try-error")) {
        
        conv <- fit.bnm_p$convergence  # 0 is success
        sauc <- SAUC(par = fit.bnm_p$par)  # u1 u2 t1 t2 r1
        bnm.p[c(1,2,4,5,7,15,16)]  <- c(fit.bnm_p$par, sauc, conv)
          
        }  
      } 
      
    ## 2 BNM.O ----
    
    {
      fn.bnm.ml  <- function(par) llk.BNM.ml(par, o.y1, o.y2, o.s11, o.s22, o.s12)
      
      fit.bnm_o  <- nlminb(bnm.int, fn.bnm.ml,
                           lower = c(rep(-Inf,2), rep(0, 2),  -1),
                           upper = c(rep( Inf,2), rep(Inf, 2), 0) 
                           )

      if(!inherits(fit.bnm_o, "try-error")) {
        
        conv  <- fit.bnm_o$convergence
        sauc  <- SAUC(par = fit.bnm_o$par)  # u1 u2 t1 t2 r1
        bnm.o[c(1,2,4,5,7,15,16)]  <- c(fit.bnm_o$par, sauc, conv)
      } 
    }
    

    ## 3. TNM.P  ----
    
    {
      fn.tnm.ml  <- function(par) llk.TNM.ml(par, p.y1, p.y2, p.y3, p.s11, p.s22, p.s33, p.s12, p.s13, p.s23)

      fit.tnm_p  <- nlminb(tnm.int, fn.tnm.ml,
                           lower = c(rep(-Inf,3), rep(0, 3),-1, rep(-1,2)),
                           upper = c(rep( Inf,3), rep(5, 3), 0, rep(1,2))
      )

      if(!inherits(fit.tnm_p, "try-error")) tnm.intp <- c(fit.tnm_p$par, 1) else tnm.intp <- c(tnm.int, 1)

    }
    
    ## 4. TNM.O  ----
    
    {
      fn.tnm.ml  <- function(par) llk.TNM.ml(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23)

      fit.tnm_o  <- nlminb(tnm.int, fn.tnm.ml,
                           lower = c(rep(-Inf,3), rep(0, 3),-1, rep(-1,2)),
                           upper = c(rep( Inf,3), rep(5, 3), 0, rep(1,2))
      )

      if(!inherits(fit.tnm_o, "try-error")) tnm.into <- c(fit.tnm_o$par, 0.5) else tnm.into <- c(tnm.int, 0.5)

    }
    
    ## 3 TNM.SA  ----
    {
      fn.tnm.ml  <- function(par) clk.TNM.ml(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p=p.hat, a.interval = c(-10, 10))

      #tnm.intp
      fit.tnm_sa_p <- nlminb(
        tnm.intp, fn.tnm.ml,
        lower = c(rep(-Inf,2), eps, rep(eps, 3), rep(-1,3), eps),
        upper = c(rep( Inf,2), Inf, rep(Inf, 3), rep( 1,3), 2) )
      
      if(!inherits(fit.tnm_sa_p, "try-error")) {
        
        conv  <- fit.tnm_sa_p$convergence
        sauc  <- SAUC(par = fit.tnm_sa_p$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
        
        ## CALCULATE alpha.p
        # b <- 1
        b   <- fit.tnm_sa_p$par[10]
        u3  <- fit.tnm_sa_p$par[3]
        t33 <- fit.tnm_sa_p$par[9]
        
        a.p <- function(a) {sum(1/f.b(a, b=b, u3=u3, t33=t33, v3=o.s33), na.rm = TRUE) - NN/p.hat}
        a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-10,10), extendInt="yes"), silent = TRUE))
        alpha.p <- a.opt.try$root
        
        p.est <- mean(pnorm(alpha.p+b*(t.p)))
        
        tnm.sa_p[c(1:11, 14:16)] <- c(fit.tnm_sa_p$par, alpha.p, p.est, sauc, conv)
          
      }  
      
      #tnm.into
      fit.tnm_sa_o <- nlminb(
        c(tnm.into), fn.tnm.ml,
        lower = c(rep(-Inf,2), eps, rep(eps, 3), rep(-1,3), eps),
        upper = c(rep( Inf,2), Inf, rep(Inf, 3), rep( 1,3), 2) )
      
      if(!inherits(fit.tnm_sa_o, "try-error")) {
        
        conv  <- fit.tnm_sa_o$convergence
        sauc  <- SAUC(par = fit.tnm_sa_o$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
        
        ## CALCULATE alpha.p
        # b <- 1
        b   <- fit.tnm_sa_o$par[10]
        u3  <- fit.tnm_sa_o$par[3]
        t33 <- fit.tnm_sa_o$par[9]
        
        a.p <- function(a) {sum(1/f.b(a, b=b, u3=u3, t33=t33, v3=o.s33), na.rm = TRUE) - NN/p.hat}
        a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-10,10), extendInt="yes"), silent = TRUE))
        alpha.p <- a.opt.try$root
        
        p.est <- mean(pnorm(alpha.p+b*(t.p)))
        
        tnm.sa_o[c(1:11, 14:16)] <- c(fit.tnm_sa_o$par, alpha.p, p.est, sauc, conv)
        
      }  
      
      #tnm.int
      fit.tnm_sa <- nlminb(
        c(tnm.int, 0.5), fn.tnm.ml,
        lower = c(rep(-Inf,2), eps, rep(eps, 3), rep(-1,3), eps),
        upper = c(rep( Inf,2), Inf, rep(Inf, 3), rep( 1,3), 2) )
      
      if(!inherits(fit.tnm_sa, "try-error")) {
        
        conv  <- fit.tnm_sa$convergence
        sauc  <- SAUC(par = fit.tnm_sa$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
        
        ## CALCULATE alpha.p
        
        b   <- fit.tnm_sa$par[10]
        u3  <- fit.tnm_sa$par[3]
        t33 <- fit.tnm_sa$par[9]
        
        a.p <- function(a) {sum(1/f.b(a, b=b, u3=u3, t33=t33, v3=o.s33), na.rm = TRUE) - NN/p.hat}
        a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-10,10), extendInt="yes"), silent = TRUE))
        alpha.p <- a.opt.try$root
        
        p.est <- mean(pnorm(alpha.p+b*(t.p)))
        
        tnm.sa[c(1:11, 14:16)] <- c(fit.tnm_sa$par, alpha.p, p.est, sauc, conv)
        
      }  
      
      
      } 
      

    }
  
  
  res <- cbind(bnm.p, bnm.o, tnm.sa, tnm.sa_o, tnm.sa_p)
  rownames(res) <- c(
    "mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3", 
    "beta", "alpha.p", "S", "N", "prop", "sauc", "conv")
  res
  
  }
  
  #######
  
  
}



