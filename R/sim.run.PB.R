##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu.run.PB <- function(S, tK, b, a, Unif.min, Unif.max, CD, Exp.r, logm, logs, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd){

  
  # files.sources <- list.files(path = "funs/"); x <- sapply(paste0("funs/", files.sources), source); library(survival); library(mixmeta)
  
  # S <- 40; tK <- 3; b <- 2; a <- -3; Unif.min <- 50; Unif.max <- 150; CD = "EXP"; Exp.r <- 0.2; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; p <- 0.7
  
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
    
    ## 2. CREATE OBSERVED LIST ----
    
    # b <- 1; a <- -2
    
    dataSt.p$slt.id <- select.id(dataSt.p, b, a)
    dataSt.o <- dataSt.p[dataSt.p$slt.id == 1, ]
    
    ## MARGINAL p
    prop <- nrow(dataSt.o)/nrow(dataSt.p)
    p.hat <- mean(pnorm(b * dataSt.p$t_lnHR + a), na.rm = TRUE)
    

    ## INITIAL b
    
    # fit1 <- glm(slt.id ~ t_lnHR, data = dataSt.p, family = binomial(link = "logit") )
    fit2 <- suppressWarnings(glm(slt.id ~ t_lnHR, data = dataSt.p, family = binomial(link = "probit") ))

    # if(fit1$converged) b0 <- unname(fit1$coefficients[2]) else b0 <- runif(1, 0.1, 1)
    if(fit2$converged) b0 <- unname(fit2$coefficients[2]) else b0 <- 0.1
    # b0 <- 0
    
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
    
    SS <- nrow(dataSt.p)
    NN <- nrow(dataSt.o)
    
    bnm.p  <- c(rep(NA, 11), SS, NA)
    tnm.p  <- c(rep(NA, 11), SS, NA)
    tnm.p2 <- c(rep(NA, 11), SS, NA)
    
    bnm.o <- c(rep(NA, 11), NN, NA)
    tnm.o <- c(rep(NA, 11), NN, NA)
    tnm.o2 <- c(rep(NA, 11), NN, NA)
    
    tnm.sa1   <- c(rep(NA, 11), prop, NA)
    tnm.sa2   <- c(rep(NA, 11), p.hat, NA)
    # tnm.sa3   <- c(rep(NA, 11), p.hat, NA)
    
    # names(bnm.par.p) <- names(tnm.p) <- 
    #   names(bnm.o) <- names(tnm.o) <- 
    #   names(tnm.sa1) <- names(tnm.sa2) <- 
    #   c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b", "sauc", "S/N/p", "conv")
    
    ## 1 BNM.P (mixmeta, bnm.p)----
    {
      Y1   <- as.matrix( cbind(p.y1, p.y2) )
      Sig1 <- as.matrix( cbind(p.s11, p.s12, p.s22) )
      
      fit.bnm_p <- try(mixmeta(Y1, Sig1, method = "ml"), silent = TRUE)
      
      if(!inherits(fit.bnm_p, "try-error")) {
        
        conv_p  <- fit.bnm_p$converged
        
        if(conv_p){
          
          u    <- fit.bnm_p$coefficients
          tt   <- c(fit.bnm_p$Psi)
          t12  <- sqrt(tt[c(1, 4)])
          r1   <- tt[2]/prod(t12)
          
          sauc  <- SAUC(par = c(u, t12, r1))  # u1 u2 t1 t2 r1
          
          bnm.p <- c(u, NA, t12, NA, r1, NA, NA, NA, sauc, SS, (1-conv_p))
          
        }  
      } 
      
      
    }
    
    ## 2 TNM.P (mixmeta tnm.p)----
  {
    Y2   <- as.matrix( cbind(p.y1, p.y2, p.y3) )
    Sig2 <- as.matrix( cbind(p.s11, p.s12, p.s13, p.s22, p.s23, p.s33) )

    fit.tnm_p <- try(mixmeta(Y2, Sig2, method = "ml"), silent = TRUE)

    if(!inherits(fit.tnm_p, "try-error")) {

      conv_p <- fit.tnm_p$converged

      if(conv_p){

        u    <- fit.tnm_p$coefficients
        tt   <- c(fit.tnm_p$Psi)
        t123 <- sqrt(tt[c(1, 5, 9)])
        r1   <- tt[2]/prod(t123[1:2])
        r2   <- tt[6]/prod(t123[2:3])
        r3   <- tt[3]/prod(t123[-2])

        sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))

        tnm.p <- c(u, t123, r1, r2, r3, NA, sauc, SS, (1-conv_p))

      }
    }

  }
    
    ## 2.1 TNM.P (expand tnm.p2)----
    {  
      int_p <- tnm.p[1:9]
      fn_p.exp <- function(par) llk.TNM.exp(par, p.y1, p.y2, p.y3, p.s11, p.s22, p.s33, p.s12, p.s13, p.s23)
      
      opt <- try(
        nlminb(int_p,
               fn_p.exp,
               lower = c(rep(-Inf,3), rep(.Machine$double.eps,3), rep(-1,3)),
               upper = c(rep( Inf,3), rep( Inf,3), rep( 1,3)) 
        ), silent = TRUE)
      
      if(!inherits(opt, "try-error")) {
        
        if(opt$convergence == 0) {
          
          
          sauc <- SAUC(par = opt$par[c(1,2,4,5,7)] )
          
          tnm.p2 <- c(opt$par, NA, sauc, SS, opt$convergence)
          
        } 
        
      } 
      
    }
    
    ## OBSERVED DATA ----
    
    
    ## 3 BNM.O (mixmeta, bnm.o) ----
    
    {
      Y3   <- as.matrix( cbind(o.y1, o.y2) )
      Sig3 <- as.matrix( cbind(o.s11, o.s12, o.s22) )
      
      fit.bnm_o <- try(mixmeta(Y3, Sig3, method = "ml"), silent = TRUE)
      
      if(!inherits(fit.bnm_o, "try-error")) {
        
        conv_o  <- fit.bnm_o$converged
        
        if(conv_o) {
          
          u    <- fit.bnm_o$coefficients
          tt   <- c(fit.bnm_o$Psi)
          t12  <- sqrt(tt[c(1, 4)])
          r1   <- tt[2]/prod(t12)
          
          sauc  <- SAUC(par = c(u, t12, r1))  # u1 u2 t1 t2 r1
          
          bnm.o <- c(u, NA,  t12, NA, r1, NA, NA, NA, sauc, NN, (1-conv_o))
          
        } 
        
      } 
      
    }
    
    ## 4 TNM.O (mixmeta, tnm.o) ----
    
    {
      Y4   <- as.matrix( cbind(o.y1, o.y2, o.y3) )
      Sig4 <- as.matrix( cbind(o.s11, o.s12, o.s13, o.s22, o.s23, o.s33) )
      
      fit.tnm_o <- try(mixmeta(Y4, Sig4, method = "ml"), silent = TRUE)
      
      if(!inherits(fit.tnm_o, "try-error")) {
        
        conv_o <- fit.tnm_o$converged
        
        if(conv_o){
          
          u    <- fit.tnm_o$coefficients
          tt   <- c(fit.tnm_o$Psi)
          t123 <- sqrt(tt[c(1, 5, 9)])
          r1   <- tt[2]/prod(t123[1:2])
          r2   <- tt[6]/prod(t123[2:3])
          r3   <- tt[3]/prod(t123[c(1,3)])
          
          sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
          
          tnm.o <- c(u, t123, r1, r2, r3, NA, sauc, NN, (1-conv_p))
          
        }
      }
      
    }
    
    ## 4.1 TNM.O (expand tnm.o2)----
    {  
      int_o <- tnm.o[1:9]
      fn_o.exp <- function(par) llk.TNM.exp(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23)
      
      opt <- try(
        nlminb(int_p,
               fn_o.exp,
               lower = c(rep(-Inf,3), rep(.Machine$double.eps,3), rep(-1,3)),
               upper = c(rep( Inf,3), rep( Inf,3), rep( 1,3)) 
        ), silent = TRUE)
      
      if(!inherits(opt, "try-error")) {
        
        if(opt$convergence == 0) {
          
          
          sauc <- SAUC(par = opt$par[c(1,2,4,5,7)] )
          
          tnm.o2 <- c(opt$par, NA, sauc, NN, opt$convergence)
          
        } 
        
      } 
      
    }
    
    
    ## SA (nlminb)----
    
    ## GET INITIAL VALUES
    # {
    #   Y4   <- as.matrix( cbind(p.y1, p.y2, p.y3) )
    #   Sig4 <- as.matrix( cbind(p.s11, p.s12, p.s13, p.s22, p.s23, p.s33) )
    #   
    #   fit.tnm_o <- try(mixmeta(Y4, Sig4, method = "ml"), silent = TRUE)
    #   
    #   if(!inherits(fit.tnm_o, "try-error")) {
    #     
    #     conv_o  <- fit.tnm_o$converged
    #     
    #     if(conv_o) {
    #       
    #       u_o    <- fit.tnm_o$coefficients
    #       tau_o  <- c(fit.tnm_o$Psi)
    #       
    #       t_o.try    <- try(c(Matrix::chol(fit.tnm_o$Psi))[c(1,4,7,5,8,9)], silent = TRUE)
    #       if (!inherits(t_o.try, "try-error")) t_o <- c(Matrix::chol(fit.tnm_o$Psi))[c(1,4,7,5,8,9)] else t_o <- rep(0.1,6)
    #       
    #       t123_o <- sqrt(tau_o[c(1, 5, 9)])
    #       r1_o   <- prod(tau_o[2]/prod(sqrt(tau_o[c(1, 5)]))) # t12
    #       r3_o   <- prod(tau_o[3]/prod(sqrt(tau_o[c(1, 9)]))) # t13
    #       r2_o   <- prod(tau_o[6]/prod(sqrt(tau_o[c(5, 9)]))) # t23
    #       
    #       
    #       init.list <- list(
    #         
    #         init01 = c(u_o, t_o),
    #         init02 = c(u_o, t123_o, r1_o, r2_o, r3_o))
    #       
    #     } else init.list <- NULL
    #     
    #   } else init.list <- NULL
    #   
    # }
    
    ## SET INITIAL VALUES
    
      # if(is.null(init.list)) int01 <- c(rep(0, 3), rep(0.1, 6), b = 2) else int01 <- c(init.list$init01, b = 2)
      # if(is.null(init.list)) int02 <- c(rep(0, 3), rep(0.1, 3), c(0.1, -0.1, 0.1), b = 2) else int02 <- c(init.list$init02, b = 2)
      # 
    if (!NA %in% tnm.o[1:9]) {
      
      # t_o.try    <- try(c(Matrix::chol(fit.tnm_o$Psi))[c(1,4,7,5,8,9)], silent = TRUE)
      # if (!inherits(t_o.try, "try-error")) t_o <- c(Matrix::chol(fit.tnm_o$Psi))[c(1,4,7,5,8,9)] else t_o <- rep(0.1,6)
      int01 <- c(round(fit.tnm_o$coefficients, 2), rep(0.1, 6))
    } else int01 <- rep(0.1, 9)
   
    if (!NA %in% tnm.o[1:9]) int02 <-  round(tnm.o[1:9], 2) else int02 <- c(rep(0.1, 3), rep(0.1, 3), c(-0.5, -0.5, 0.1))
    # int0 <- rep(0, 10)
    
    ## 4 TNM.SA1 (cholesky, tnm.sa1) ----
    {  
    fn.cho <- function(par) llk.TNMO.cho(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p = p.hat, a.interval= c(-10, 10))
    
    opt.sa1 <- try(
      nlminb(c(int01, b=0.5),
             fn.cho,
             lower = c(rep(-Inf,9), 0),
             upper = c(rep( Inf,9), 2) 
             ), silent = TRUE)
    
    if(!inherits(opt.sa1, "try-error")) {
      
      if(opt.sa1$convergence == 0) {
        
        u <- opt.sa1$par[1:3]
        
        L <- diag(rep(0,3))
        L[lower.tri(L, diag = TRUE)] <- opt.sa1$par[c(4:9)]
        Omega <- L %*% t(L)
        
        tt   <- c(Omega)
        t123 <- sqrt(tt[c(1, 5, 9)])
        r1   <- tt[2]/prod(t123[1:2])
        r2   <- tt[6]/prod(t123[2:3])
        r3   <- tt[3]/prod(t123[c(1,3)])
        
        sauc <- SAUC(par = c(u[1:2], t123[1:2], r1))
        b    <- opt.sa1$par[10]
        tnm.sa1 <- c(u, t123, r1, r2, r3, b, sauc, p.hat, opt.sa1$convergence)
        
        
        # f.b <- function(a){
        #   
        #   sq <- suppressWarnings(sqrt(1 + b^2 * (1 + Omega[3,3] / o.s33)))
        #   
        #   pnorm( (a + b * u[3] / sqrt(o.s33)) / sq )
        #   
        # }
        # 
        # 
        # a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p.hat}
        # 
        # a.opt.try <- suppressWarnings(try(uniroot(a.p, interval = c(-10,10), extendInt="yes"), silent = FALSE))
        # 
        # a.opt <- a.opt.try$root
        
        
      } 
      
    } 
    

    }
      ## 5. TNM.SA2 (expand, tnm.sa2) ----
      {  
        fn.exp <- function(par) llk.TNMO.exp(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p = p.hat, a.interval= c(-10, 10))
        
        opt.sa2 <- try(
          nlminb(c(int02, b=0.5),
                 fn.exp,
                 lower = c(rep(-Inf,3), rep(.Machine$double.eps,3), rep(-1,3), 0),
                 upper = c(rep( Inf,3), rep( 3,3), rep( 1,3), 2) 
          ), silent = TRUE)
        
        if(!inherits(opt.sa2, "try-error")) {
          
          if(opt.sa2$convergence == 0) {

            
            sauc <- SAUC(par = opt.sa2$par[c(1,2,4,5,7)] )
            
            tnm.sa2 <- c(opt.sa2$par, sauc, prop, opt.sa2$convergence)
            
          } 
          
        } 

      }
      
    
    # {  
    #   fn.exp <- function(par)llk.TNM.o(par, o.y1, o.y2, o.y3, o.s11, o.s22, o.s33, o.s12, o.s13, o.s23, p = p.hat, a.interval= c(-10, 10))
    #   
    #   opt.sa2 <- try(
    #     nlminb(c(int02, b=0.5),
    #            fn.exp,
    #            lower = c(rep(-Inf,3), rep(.Machine$double.eps,3), rep(-1,3), 0),
    #            upper = c(rep( Inf,3), rep( 3,3), rep( 1,3), 2) 
    #     ), silent = TRUE)
    #   
    #   if(!inherits(opt.sa2, "try-error")) {
    #     
    #     if(opt.sa2$convergence == 0) {
    #       
    #       
    #       sauc <- SAUC(par = opt.sa2$par[c(1,2,4,5,7)] )
    #       
    #       tnm.sa3 <- c(opt.sa2$par, sauc, prop, opt.sa2$convergence)
    #       
    #     } 
    #     
    #   } 
    #   
    # }
  
  res <- cbind(bnm.p, tnm.p, tnm.p2, bnm.o, tnm.o, tnm.o2, tnm.sa1, tnm.sa2)
  rownames(res) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b", "sauc", "S/N/p", "conv")
  res
  
  }
  
  
  #######
  
  
}



