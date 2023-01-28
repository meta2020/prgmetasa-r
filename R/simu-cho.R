##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu.cho <- function(S, tK, b, a, N.min, N.max, CD, Exp.r, meanlog, sdlog, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd){

  
  # files.sources <- list.files(path = ""); x <- sapply(paste0("", files.sources), source); 
  # library(survival); library(mixmeta)
  # load(sprintf("a-calc/beta_1_alpha_%d.RData", N))
  # S <- 40; tK <- 3; b <- 2; a <- -3; N.min <- 50; N.max <- 150; CD = "EXP"; Exp.r <- 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; p <- 0.7
  # S = N/0.8; tK = 2; b = 1.5; a = set$alpha1; N.min = 10; N.max= 50; CD = "EXP"; Exp.r = 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2

  eps <- sqrt(.Machine$double.eps)
  
  ####### RUN
  {
    ## POPULATION DATA
    
    p_dataSt.all <- population.data.list(
      S, tK, N.min, N.max, 
      CD, Exp.r, meanlog, sdlog, minc, maxc, 
      x1.u, x2.u, x1.s, x2.s, v.sd)
    p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))
    
    ## REMOVE NA
    
    p_dataSt <- p_dataSt[complete.cases(p_dataSt),]
    
    ## ANALYTICAL POPULATION DATA 
    
    eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
    dataSt.p <- convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)
    
    ## REMOVE NA
    
    dataSt.p <- dataSt.p[complete.cases(dataSt.p),]
    
    ## 2. CREATE OBSERVED LIST ----
    
    # b <- 3; a <- -1
    
    dataSt.p$slt.id <- select.id(dataSt.p, b, a)
    dataSt.o <- dataSt.p[dataSt.p$slt.id == 1, ]
    
    
    SS <- nrow(dataSt.p)
    NN <- nrow(dataSt.o)
    
    ## MARGINAL p
    prop <- nrow(dataSt.o)/nrow(dataSt.p)
    p.hat <- mean(pnorm(b * dataSt.p$t_lnHR + a), na.rm = TRUE)
    
    bnm.p  <- c(rep(NA, 11), SS, NN, prop, NA, NA)
    bnm.o <- tnm.p <- tnm.o <-tnm.sa <- rep(NA, 16)
    names(bnm.p) <- names(bnm.o) <- names(tnm.p) <- names(tnm.o) <- names(tnm.sa) <- c(
      "mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3", 
      "beta", "alpha", "S", "N", "prop", "sauc","conv")
    
    
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

    bnm.int <- c(rep(0.1,4), -0.5)
    tnm.int <- rep(0.1,9)
    b0 <- 0.5
    
    
    ## 1 BNM.P ----
    {
      
      fn.bnm.ml  <- function(par) llk.BNM.ml(par, p.y1, p.y2, p.s11, p.s22, p.s12)
      fit.bnm_p  <- nlminb(bnm.int, fn.bnm.ml,
                           lower = c(rep(-Inf,2), rep(0, 2),  -1),
                           upper = c(rep( Inf,2), rep(Inf, 2), 1) 
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
                           upper = c(rep( Inf,2), rep(Inf, 2), 1) 
                           )

      if(!inherits(fit.bnm_o, "try-error")) {
        
        conv  <- fit.bnm_o$convergence
        sauc  <- SAUC(par = fit.bnm_o$par)  # u1 u2 t1 t2 r1
        bnm.o[c(1,2,4,5,7,15,16)]  <- c(fit.bnm_o$par, sauc, conv)
      } 
    }
    

    ## 3. TNM.P  ----
    
    {
      fn.tnm.ml  <- function(par) llk.TNM.cho(par, p.y1, p.y2, p.y3, p.s11, p.s22, p.s33, p.s12, p.s13, p.s23)
      
      fit.tnm_p  <- nlminb(tnm.int, fn.tnm.ml,
                           lower = c(rep(-Inf,3), rep(eps, 6)),
                           upper = c(rep( Inf,3), rep(Inf, 6)) 
      )
      
      if(!inherits(fit.tnm_p, "try-error")) {
        
        conv  <- fit.tnm_p$convergence
        mu <- fit.tnm_p$par[1:3]
        L <- diag(rep(0,3))
        L[lower.tri(L, diag = TRUE)] <- fit.tnm_p$par[4:9]
        O <- L %*% t(L)
        tt <- c(O)
        t123 <- sqrt(tt[c(1,5,9)])
        r1 <- tt[2]/prod(t123[c(1,2)])
        r2 <- tt[6]/prod(t123[c(2,3)])
        r3 <- tt[3]/prod(t123[c(1,3)])
        
        sauc  <- SAUC(par = c(mu[1:2], t123[1:2], r1))  # u1 u2 t1 t2 r1
        tnm.p[c(1:9, 15,16)]  <- c(mu, t123, r1, r2, r3, sauc, conv)
        
      }
      
    }
    
    ## 4. TNM.O  ----
    
    {
      fn.tnm.ml  <- function(par) llk.TNM.cho(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23)
      
      fit.tnm_o  <- nlminb(tnm.int, fn.tnm.ml,
                           lower = c(rep(-Inf,3), rep(0, 3),  rep(-1,3)),
                           upper = c(rep( Inf,3), rep(5, 3), rep(1,3)) 
      )
      
      if(!inherits(fit.tnm_o, "try-error")) {
        
        conv  <- fit.tnm_o$convergence
        mu <- fit.tnm_o$par[1:3]
        L <- diag(rep(0,3))
        L[lower.tri(L, diag = TRUE)] <- fit.tnm_o$par[4:9]
        O <- L %*% t(L)
        tt <- c(O)
        t123 <- sqrt(tt[c(1,5,9)])
        r1 <- tt[2]/prod(t123[c(1,2)])
        r2 <- tt[6]/prod(t123[c(2,3)])
        r3 <- tt[3]/prod(t123[c(1,3)])
        
        sauc  <- SAUC(par = c(mu[1:2], t123[1:2], r1))  # u1 u2 t1 t2 r1
        tnm.o[c(1:9, 15,16)]  <- c(mu, t123, r1, r2, r3, sauc, conv)
        
      }
      
    }
    ## 3 TNM.SA  ----
    
    {
      fn.tnm.ml  <- function(par) clk.TNM.cho(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p=p.hat, a.interval = c(-5, 5))
      
      fit.tnm_sa <- nlminb(
        c(tnm.int, b0), fn.tnm.ml,
        lower = c(rep(-Inf,3), rep(eps, 3), rep(-1,3), eps),
        upper = c(rep( Inf,3), rep(Inf, 3), rep( 1,3), 2) )
      
      if(!inherits(fit.tnm_sa, "try-error")) {
        
        conv  <- fit.tnm_sa$convergence
        mu <- fit.tnm_sa$par[1:3]
        L <- diag(rep(0,3))
        L[lower.tri(L, diag = TRUE)] <- fit.tnm_sa$par[4:9]
        O <- L %*% t(L)
        tt <- c(O)
        t123 <- sqrt(tt[c(1,5,9)])
        r1 <- tt[2]/prod(t123[c(1,2)])
        r2 <- tt[6]/prod(t123[c(2,3)])
        r3 <- tt[3]/prod(t123[c(1,3)])
        
        sauc  <- SAUC(par = c(mu[1:2], t123[1:2], r1))  # u1 u2 t1 t2 r1
        
        ## CALCULATE ALPHA
        b   <- fit.tnm_sa$par[10]
        u3  <- mu[3]
        t33 <- t123[3]^2
        v3  <- o.s33
        
        
        f.b <- function(a){
          
          sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / v3)))
          
          pnorm( (a + b * u3/sqrt(v3)) / sq )
          
        }
        a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p.hat}
        
        a.opt.try <- suppressWarnings(try(uniroot(a.p, c(-5,5), extendInt="no"), silent = TRUE)) 
        
        alpha <- a.opt.try$root
        
        tnm.sa[c(1:11, 15,16)] <- c(mu, t123, r1, r2, r3, b, alpha, sauc, conv)
          
        }  
      } 
      
      
    }
  
  
  res <- cbind(bnm.p, bnm.o, tnm.p, tnm.o, tnm.sa)
  res
  
  }
  
  #######
  
  
}



