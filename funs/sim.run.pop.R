##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu.run.pop <- function(S, tK, b, a, Unif.min, Unif.max, CD, Exp.r, logm, logs, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd){

  
  # files.sources <- list.files(path = "funs/"); x <- sapply(paste0("funs/", files.sources), source); library(survival);library(mixmeta)
  # S <- 30; tK <- 3; b <- 0.5; a <- -0.7; Unif.min <- 50; Unif.max <- 100; CD = "EXP"; Exp.r <- 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; p <- 0.7

  ####### RUN
  {
    ## POPULATION DATA
    
    p_dataSt.all <- population.data.list(S, tK, Unif.min, Unif.max, CD, Exp.r, logm, logs, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd)
    p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))
    
    ## REMOVE NA
    
    p_dataSt <- p_dataSt[complete.cases(p_dataSt),]
    
    ## ANALYTICAL POPULATION DATA 
    eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
    dataSt.p <- convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)
    
    ## REMOVE NA
    
    dataSt.p <- dataSt.p[complete.cases(dataSt.p),]
    
    
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

    }

    ## 3. OPTIMIZE POPULATION----
    
   
    tnm.mix <- tnm.cho <- tnm.exp <- bnm.mix <- bnm.exp <- rep(NA, 11) 
    # names(tnm.mix) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3")
    
    ## 3.1 TNM.P (mixmeta)----
  {    
    Y1   <- as.matrix( cbind(p.y1, p.y2, p.y3) )
    Sig1 <- as.matrix( cbind(p.s11, p.s12, p.s13, p.s22, p.s23, p.s33) )
    
    fit.tnm_p <- try(mixmeta(Y1, Sig1, method = "ml"), silent = TRUE)
    
    if(!inherits(fit.tnm_p, "try-error")) {
      
      conv_p <- fit.tnm_p$converged
      
      if(conv_p){
        
        u    <- fit.tnm_p$coefficients
        tt   <- c(fit.tnm_p$Psi)
        t123 <- sqrt(tt[c(1, 5, 9)])  # t1 t2 t3
        r1   <- tt[2]/prod(t123[1:2]) # t12/(t1*t2)
        r2   <- tt[6]/prod(t123[2:3])
        r3   <- tt[3]/prod(t123[-2])
        
        sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
        
        tnm.mix <- c(u, t123, r1, r2, r3, sauc, (1-conv_p))
        
      } 
      
    }
    
  }
   
    ## 3.2 TNM.P (TNM.cho)----
    {    

      fn.tnm <- function(par) llk.TNM.cho(par, y1 = p.y1, y2 = p.y2, y3 = p.y3, s11 = p.s11, s22 = p.s22, s33 = p.s33, s12 = p.s12, s13 = p.s13, s23 = p.s23)
      int0 <- c(rep(0, 3), rep(0.1, 6))
      
      # opt <- try(nlminb(int0,
      #               fn.tnm,
      #               lower = c(rep(-Inf,3), rep(-Inf,6)),
      #               upper = c(rep( Inf,3), rep(Inf,6))
      # ), silent = TRUE)
      
      opt <- try(optim(int0,
                        fn.tnm,
                       method = "L-BFGS-B"
      ), silent = TRUE)
      
      if(!inherits(opt, "try-error")) {
        
        conv_p <- opt$convergence
        
        if(conv_p == 0){
          
          L <- diag(rep(0,3))
          L[lower.tri(L, diag = TRUE)] <- opt$par[4:9]
          Omega <- L %*% t(L)
          
          u    <- opt$par[1:3]
          tt   <- c(Omega)
          t123 <- sqrt(tt[c(1, 5, 9)])  # t1 t2 t3
          r1   <- tt[2]/prod(t123[1:2]) # t12/(t1*t2)
          r2   <- tt[6]/prod(t123[2:3])
          r3   <- tt[3]/prod(t123[-2])
          
          sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
          
          tnm.cho <- c(u, t123, r1, r2, r3, sauc, conv_p)
          
        } 
        
      } 
      
    } 
   

    # ## 3.3 TNM.P (TNM.matrix, NO NEED)----
    # {    
    #   
    #   fn.tnm <- function(par) llk.TNM.mat(par, y1 = p.y1, y2 = p.y2, y3 = p.y3, s11 = p.s11, s22 = p.s22, s33 = p.s33, s12 = p.s12, s13 = p.s13, s23 = p.s23)
    #   int0 <- c(rep(0, 3), rep(0.1, 6))
    #   
    #   opt <- nlminb(int0,
    #                 fn.tnm,
    #                 lower = c(rep(-Inf,3), rep(.Machine$double.eps,3), rep(-1, 3)),
    #                 upper = c(rep( Inf,3), rep(Inf,3), rep(1,3))
    #   )
    #   
    #   if(!inherits(opt, "try-error")) {
    #     
    #     conv_p <- opt$convergence
    #     
    #     if(conv_p == 0){
    #       
    #       u_p    <- opt$par[1:3]
    #       t123_p <- (opt$par[4:6])^2
    #       r      <- opt$par[7:9]
    #       t12    <- sqrt(prod(t123_p[c(1,3)]))*r[1]
    #       t23    <- sqrt(prod(t123_p[c(2,3)]))*r[2]
    #       t13    <- sqrt(prod(t123_p[c(1,3)]))*r[3]
    #       
    #       sauc_p <- SAUC(par = opt$par[c(1,2,4,5,7)])
    #       
    #       tnm.par.p.mat <- c(u_p, t123_p, t12, t23, t13, sauc_p, conv_p)
    #       
    #     } else tnm.par.p.mat <- rep(NA, 11) 
    #     
    #   } else tnm.par.p.mat <- rep(NA, 11)
    #   
    #   names(tnm.par.p.mat) <- c("u1", "u2", "u3", "t11", "t22", "t33", "t12", "t23", "t13","sauc", "conv")
    # } 
    # 
    
    ## 3.4 TNM.P (TNM.expand)----
    {    
      
      fn.tnm <- function(par) llk.TNM.exp(par, y1 = p.y1, y2 = p.y2, y3 = p.y3, s11 = p.s11, s22 = p.s22, s33 = p.s33, s12 = p.s12, s13 = p.s13, s23 = p.s23)
      int0 <- c(rep(0, 3), rep(0.1, 3), c(-0.1, -0.1, 0.1))
      
      opt <- try(nlminb(int0,
                    fn.tnm,
                    lower = c(rep(-5,3), rep(.Machine$double.eps,3), rep(-1, 3)),
                    upper = c(rep( 5,3), rep(3,3), rep(1,3))
      ), silent = TRUE)
      
      if(!inherits(opt, "try-error")) {
        
        conv_p <- opt$convergence
        
        if(conv_p == 0){

          sauc_p <- SAUC(par = opt$par[c(1,2,4,5,7)])
          
          tnm.exp <- c(opt$par, sauc_p, conv_p)
          
        } 
        
      } 
      
    } 
    
    ## 3.5 BNM.P (mixmeta)----
    {    
      Y5   <- as.matrix( cbind(p.y1, p.y2) )
      Sig5 <- as.matrix( cbind(p.s11, p.s12, p.s22) )
      
      fit.bnm_p <- try(mixmeta(Y5, Sig5, method = "ml"), silent = TRUE)
      
      if(!inherits(fit.bnm_p, "try-error")) {
        
        conv_p <- fit.bnm_p$converged
        
        if(conv_p){
          
          u    <- c(fit.bnm_p$coefficients, NA)
          tt   <- c(fit.bnm_p$Psi)
          t123 <- c(sqrt(tt[c(1, 4)]), NA)  # t1 t2 t3
          r1   <- tt[2]/prod(t123[1:2]) # t12/(t1*t2)
          r2   <- NA
          r3   <- NA
          
          sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
          
          bnm.mix <- c(u, t123, r1, r2, r3, sauc, (1-conv_p))
          
        } 
        
      }
      
    }
    
    
    ## 3.6 BNM.P (BNM.expand)----
    {    
      
      fn.bnm <- function(par) llk.BNM.exp(par, y1 = p.y1, y2 = p.y2, v1 = p.s11, v2 = p.s22, v12 = p.s12)
      int0 <- c(rep(0, 2), rep(0.1, 2), c(-0.1))
      
      opt <- try(nlminb(int0,
                        fn.bnm,
                        lower = c(rep(-Inf,2), rep(.Machine$double.eps,2), -1),
                        upper = c(rep( Inf,2), rep(3,2), 1)
      ), silent = TRUE)
      
      if(!inherits(opt, "try-error")) {
        
        conv_p <- opt$convergence
        
        if(conv_p == 0){
          
          u    <- c(opt$par[1:2], NA)
          t123 <- c(opt$par[3:4], NA)
          r1   <- opt$par[5]
          r2   <- NA
          r3   <- NA
          
          sauc <- SAUC(par = opt$par)
          
          bnm.exp <- c(u, t123, r1, r2, r3, sauc, conv_p)
          
        } 
        
      } 
      
    } 
    
    
    
  ## RESULTS
    
  res <- cbind(tnm.mix, 
               tnm.cho, 
               # tnm.par.p.mat, 
               tnm.exp,
               bnm.mix, 
               bnm.exp)
  rownames(res) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "sauc", "conv")
  res
  
  }
  
  
  #######
  
  
}



