##
## SIMULATION ALL 
##


## 1. CREATE POPULATION LIST ----

simu <- function(
    S, tK, b, a, N.min, N.max, CD, Exp.r, 
    meanlog, sdlog, minc, maxc, 
    x1.u, x2.u, x1.s, x2.s, v.sd, 
    true.p=NULL){


  # {S = 200; tK = 2; b = 1; a = set$alpha2; N.min = set$pts.dist[1]; N.max= set$pts.dist[2]
  # CD = "EXP"; Exp.r = 0.2
  # x1.u = set$mk2[1]; x2.u = set$mk2[2]; x1.s = set$mk2[3]; x2.s = set$mk2[4]; v.sd = set$mk2[5]
  # true.p = NULL
  # eps <- sqrt(.Machine$double.eps)
  # }
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
    S.p <- nrow(dataSt.p)
    
    ## 2. CREATE OBSERVED LIST FOR 0=0.7, P=0.5, P=0.3 ----
    
    # P=0.7
    select.p.7  <- pnorm(b * (dataSt.p$t_lnHR) + a[1])
    select.id.7 <- rbinom(S.p, 1, select.p.7) 
    dataSt.o.7 <- dataSt.p[select.id.7 == 1, ]
    N.7 <- nrow(dataSt.o.7)
    pp.7 <- N.7/S.p
    pe.7 <- mean(select.p.7)
    if (is.null(true.p)) p.hat.7 <- pp.7  else p.hat.7 <- true.p 
    # P=0.5
    select.p.5  <- pnorm(b * (dataSt.p$t_lnHR) + a[2])
    select.id.5 <- rbinom(S.p, 1, select.p.5 ) 
    dataSt.o.5 <- dataSt.p[select.id.5 == 1, ]
    N.5 <- nrow(dataSt.o.5)
    pp.5 <- N.5/S.p
    pe.5 <- mean(select.p.5)
    if (is.null(true.p)) p.hat.5 <- N.5/S.p else p.hat.5 <- true.p 
    # P=0.3
    select.p.3  <- pnorm(b * (dataSt.p$t_lnHR) + a[3])
    select.id.3 <- rbinom(S.p, 1, select.p.3 ) 
    dataSt.o.3 <- dataSt.p[select.id.3 == 1, ]
    N.3 <- nrow(dataSt.o.3)
    pp.3 <- N.3/S.p
    pe.3 <- mean(select.p.3)
    if (is.null(true.p)) p.hat.3 <- N.3/S.p else p.hat.3 <- true.p 
    
    ## MARGINAL p

    
    # c("mu1", "mu2", "mu3", "tau1", "tau2", "tau3", "rho1", "rho2", "rho3", "sauc", "conv",
    #   "beta", "alpha.p", "SN", "p.est", "p.est.a", "p.est.b")
    
    # VECTORS TO STORE RESULTS
    bnm.p <- c(rep(NA, 13), S.p, NA, NA, NA)
    tnm.p <- rep(NA, 17)
    
    bnm.o <- rep(NA, 17)
    tnm.o <- rep(NA, 17)
    tnm.sa_p <- tnm.sa_o <-tnm.sa_n <- rep(NA, 17)
    

    ## SET INTITIAL VALUES

    # tnm.int <- c(runif(2, -1, 1), runif(4, 0.1, 1), runif(3, -1, -0.1))
    # tnm.int <- c(rep(0.1,6), rep(-0.1, 3))
    
    
    # tnm.intp <- tnm.int
    
    b0 <- 0.5
    
    
    ## 1.1 BNM.P ----
    {
      
      fn.bnm.ml  <- function(par) llk.BNM.ml(par, dataSt.p$u_sen, dataSt.p$u_spe, dataSt.p$v_sen, dataSt.p$v_spe, dataSt.p$v_senspe)
      fit.bnm_p  <- nlminb(
        c(rep(0.1,4),-0.1),  ## initial values for mu1,2 tau1,2 rho1
        fn.bnm.ml,
        lower = c(rep(-5,2), rep(0, 2),-1),
        upper = c(rep( 5,2), rep(3, 2), 1)
        )

      if(!inherits(fit.bnm_p, "try-error")) {

        conv <- fit.bnm_p$convergence  # 0 is success
        sauc <- SAUC(par = fit.bnm_p$par)  # u1 u2 t1 t2 r1
        bnm.p[c(1,2,4,5,7,10,11)]  <- c(fit.bnm_p$par, sauc, conv)

        }
      
      # MIXMETA USED CHOLESKY DECOMPOSITION
      # fit.bnm_p <- try(mixmeta::mixmeta(cbind(dataSt.p$u_sen, dataSt.p$u_spe),S=cbind(dataSt.p$v_sen, dataSt.p$v_senspe, dataSt.p$v_spe), data=dataSt.p, method = "ml"), silent = TRUE)
      # 
      # if((!inherits(fit.bnm_p, "try-error")) ) {
      #   
      #   conv <- fit.bnm_p$converged  # 0 is success
      #   mu <- fit.bnm_p$coefficients
      #   v  <- fit.bnm_p$Psi
      #   tau_rho <- c(sqrt(v[c(1,4)]), v[3]/prod(sqrt(v[c(1,4)])))
      #   sauc <- SAUC(par = c(mu, tau_rho))  # u1 u2 t1 t2 r1
      #   bnm.p[c(1,2,4,5,7,10,11)]  <- c(mu, tau_rho, sauc, conv-1)
      #   
      # } 
    } 
    
    ## 1.2 TNM.P  ----
    
    {
      fn.tnm.ml  <- function(par) llk.TNM.ml(
        par, 
        dataSt.p$u_sen, dataSt.p$u_spe, dataSt.p$u_lnHR,
        dataSt.p$v_sen, dataSt.p$v_spe, dataSt.p$v_lnHR,
        dataSt.p$v_senspe, dataSt.p$v_senlnHR, dataSt.p$v_spelnHR)
      
      # initial values
      if(!inherits(fit.bnm_p, "try-error")) {
        
        if (fit.bnm_p$convergence==0) {
          
          tnm.intp <- bnm.p[1:9]
          tnm.intp[c(3,6)] <- c(0.1,0.1)#runif(2, 0.1, 1)
          tnm.intp[c(8,9)] <- c(-0.1,-0.1)#runif(2, -1, 1)
        
      } else tnm.intp <- c(rep(0.1,6), rep(-0.1,3))} else tnm.intp <- c(rep(0.1,6), rep(-0.1,3))
      
      
      fit.tnm_p  <- nlminb(
        tnm.intp, fn.tnm.ml,
        lower = c(rep(-5,3), rep(0, 3),-1, rep(-1,2)),
        upper = c(rep( 5,3), rep(3, 3), 1, rep(1,2))
        )
      
      if(!inherits(fit.tnm_p, "try-error")) {

        conv  <- fit.tnm_p$convergence
        sauc  <- SAUC(par = fit.tnm_p$par[c(1,2,4,5,7)])  # u1 u2 t1 t2 r1
        tnm.p[1:11] <- c(fit.tnm_p$par, sauc, conv)
        
        if (conv==0) tnmsa.intp <- c(fit.tnm_p$par, 1) else tnmsa.intp <- c(tnm.intp, 1)
        
      } else tnmsa.intp <- c(tnm.intp, 1)
      
    }
    
    if(nrow(dataSt.o.7)==0) res.7 <- NULL else res.7 <- estimates.o(
      dataSt.o = dataSt.o.7, N = N.7, p.hat =p.hat.7, pp=pp.7, pe=pe.7, 
      # bnm.intoo= bnm.int, tnm.intoo=tnm.into, 
      tnmsa.intpp = tnmsa.intp, 
      bnm.oo = bnm.o, tnm.oo = tnm.o, 
      # tnm.intt = tnm.int, 
      dataSt.pp = dataSt.p)
    if(nrow(dataSt.o.5)==0) res.5 <- NULL else res.5 <- estimates.o(
      dataSt.o = dataSt.o.5, N = N.5, p.hat =p.hat.5, pp=pp.5, pe=pe.5, 
      # bnm.intoo= bnm.int, tnm.intoo=tnm.into, 
      tnmsa.intpp = tnmsa.intp, 
      bnm.oo = bnm.o, tnm.oo = tnm.o, 
      # tnm.intt = tnm.int, 
      dataSt.pp = dataSt.p)
    if(nrow(dataSt.o.3)==0) res.3 <- NULL else res.3 <- estimates.o(
      dataSt.o = dataSt.o.3, N = N.3, p.hat =p.hat.3, pp=pp.3, pe=pe.3, 
      # bnm.intoo= bnm.int, tnm.intoo=tnm.into, 
      tnmsa.intpp = tnmsa.intp, 
      bnm.oo = bnm.o, tnm.oo = tnm.o, 
      # tnm.intt = tnm.int, 
      dataSt.pp = dataSt.p)
 
  res <- cbind(bnm.p, tnm.p, res.7, res.5, res.3)
  colnames(res) <- c("BNM.P", "TNM.P", paste0(rep(c("BNM.O", "TNM.O", "SA.o", "SA.p", "SA.n"),3),rep(c(7,5,3), each=5)))
  # return(res)
  res
  }
  
     #######
  
  
}



