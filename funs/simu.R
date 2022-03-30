##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu <- function(S, tK, b, a, N.min, N.max, CD, Exp.r, logm, logs, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd){

  
  # files.sources <- list.files(path = ""); x <- sapply(paste0("", files.sources), source); library(survival); library(mixmeta)
  # load(sprintf("a-calc/beta_1_alpha_%d.RData", N))
  # S <- 40; tK <- 3; b <- 2; a <- -3; N.min <- 50; N.max <- 150; CD = "EXP"; Exp.r <- 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; p <- 0.7
  # S = N/0.2; tK = 2; b = 1; a = a[12]; N.min = 50; N.max= 150; CD = "EXP"; Exp.r = 0.2; x1.u = x1.u; x2.u = x2.u; x1.s = x1.s; x2.s = x2.s; v.sd = v.sd
  
  ####### RUN
  {
    ## POPULATION DATA
    
    p_dataSt.all <- population.data.list(S, tK, N.min, N.max, CD, Exp.r, logm, logs, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd)
    p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))
    
    ## REMOVE NA
    
    p_dataSt <- p_dataSt[complete.cases(p_dataSt),]
    
    ## ANALYTICAL POPULATION DATA 
    eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
    dataSt.p <- convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)
    
    ## REMOVE NA
    
    dataSt.p <- dataSt.p[complete.cases(dataSt.p),]
    
    ## 2. CREATE OBSERVED LIST ----
    
    # b <- 1; a <- -2
    
    dataSt.p$slt.id <- select.id(dataSt.p, b, a)
    dataSt.o <- dataSt.p[dataSt.p$slt.id == 1, ]
    
    
    SS <- nrow(dataSt.p)
    NN <- nrow(dataSt.o)
    
    ## MARGINAL p
    prop <- nrow(dataSt.o)/nrow(dataSt.p)
    p.hat <- mean(pnorm(b * dataSt.p$t_lnHR + a), na.rm = TRUE)
    
    
    bnm.p  <- c(rep(NA, 11), SS, NA)
    bnm.o  <- c(rep(NA, 11), NN, NA)
    tnm.sa <- c(rep(NA, 11), prop, NA)
    
    
    if(nrow(dataSt.o)==0) {
      
      res <- matrix(NA, nrow = 13, ncol = 3)
      rownames(res) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b", "sauc", "S/N/p", "conv")
      res
      
    } else {
      
    ## INITIAL b
    
    b0 <- 0.5
    
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

    ## POPULATION----
    
    
    tnm.int<- colMeans(cbind(o.y1, o.y2,o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23))
    
    ## 1 BNM.P ----
    {
      bnm.int    <- c(mean(p.y1), mean(p.y2), mean(p.s11), mean(p.s22), mean(p.s12))
      fn.bnm.ml  <- function(par) .llk.BNM.ml(par, p.y1, p.y2, p.s11, p.s22, p.s12)
      fit.bnm_p  <- nlminb(bnm.int, fn.bnm.ml,
                           lower = c(rep(-Inf,2), rep(0, 2),   rep(-1,2)),
                           upper = c(rep( Inf,2), rep(Inf, 2), rep( 1,2)) )
      # par.bnm.ml <- opt.bnm.ml$par

      if(!inherits(fit.bnm_p, "try-error")) {
        
        conv  <- fit.bnm_p$convergence
        
        if(conv==0){
          
          sauc  <- SAUC(par = fit.bnm_p$par)  # u1 u2 t1 t2 r1
          
          bnm.p <- c(fit.bnm_p$par[1:2], NA, fit.bnm_p$par[3:4], NA, fit.bnm_p$par[5], NA, NA, NA, sauc, SS, conv)
          
        }  
      } 
      
      
    }
    
    ## 2 BNM.O ----
    
    {
      fn.bnm.ml  <- function(par) .llk.BNM.ml(par, o.y1, o.y2, o.s11, o.s22, o.s12)
      fit.bnm_o  <- nlminb(fit.bnm_p$par, fn.bnm.ml,
                           lower = c(rep(-Inf,2), rep(0, 2),   rep(-1,2)),
                           upper = c(rep( Inf,2), rep(Inf, 2), rep( 1,2)) 
                           )
      # par.bnm.ml <- opt.bnm.ml$par
      
      if(!inherits(fit.bnm_o, "try-error")) {
        
        conv  <- fit.bnm_o$convergence
        
        if(conv==0){
          
          sauc  <- SAUC(par = fit.bnm_o$par)  # u1 u2 t1 t2 r1
          
          bnm.o <- c(fit.bnm_o$par[1:2], NA, fit.bnm_o$par[3:4], NA, fit.bnm_o$par[5], NA, NA, NA, sauc, NN, conv)
          
        }  
      } 
    }
    
    ## OBSERVED DATA ----
    
    ## Initia TNM.O (mixmeta, tnm.o) ----
    
    {
      Y   <- as.matrix( cbind(o.y1, o.y2, o.y3) )
      Sig <- as.matrix( cbind(o.s11, o.s12, o.s13, o.s22, o.s23, o.s33) )
      
      fit.tnm_o <- try(mixmeta(Y, Sig, method = "ml"), silent = TRUE)
      
      if(!inherits(fit.tnm_o, "try-error")) {
        
        if(fit.tnm_o$converged){ 
          
          u    <- fit.tnm_o$coefficients
          tt   <- c(fit.tnm_o$Psi)
          t123 <- sqrt(tt[c(1, 5, 9)])
          r1   <- tt[2]/prod(t123[1:2])
          r2   <- tt[6]/prod(t123[2:3])
          r3   <- tt[3]/prod(t123[c(1,3)])
          
          #sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
          
          tnm.int <- c(mu = unname(u), tau = t123, rho1 = r1, rho2 = r2, rho3 = r3)}
          
      }
      
    }
    
    ## 3 TNM.SA  ----
    
    {
      fn.tnm.ml  <- function(par) .clk.TNM.ml(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p=p.hat, a.interval = c(-10, 10))
      fit.tnm_sa <- nlminb(
        c(tnm.int, b0), fn.tnm.ml,
        lower = c(rep(-Inf,3), rep(0, 3),   rep(-1,3), 0),
        upper = c(rep( Inf,3), rep(Inf, 3), rep( 1,3), 2) )
      
      if(!inherits(fit.tnm_sa, "try-error")) {
        
        conv  <- fit.tnm_sa$convergence
        
        if(conv==0){
          
          sauc  <- SAUC(par = fit.tnm_sa$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
          
          tnm.sa <- c(fit.tnm_sa$par, sauc, prop, conv)
          
        }  
      } 
      
      
    }
  
  
  res <- cbind(bnm.p, bnm.o, tnm.sa)
  rownames(res) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b", "sauc", "S/N/p", "conv")
  res
  
  }
  
  }
  
  #######
  
  
}



