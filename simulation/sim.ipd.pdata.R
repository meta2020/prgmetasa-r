.SAUC <- function(par){
  
  u1  <- par[1]
  u2  <- par[2]
  t1  <- par[3]
  t2  <- par[4]
  r   <- par[5]
  
  if (NA %in% par) {auc <- NA} else {
    
    auc.try <- try(integrate(function(x) { plogis(u1 - (t1*r/t2) * (qlogis(x) + u2)) }, 0, 1))
    
    if(!inherits(auc.try,"try-error")) auc.try$value else NA
    
  }
  
}


## GENERATE IPD DATA
.p.sim.ipd <- function(
    S, 
    tK, 
    npt.range,  ## NUMBER OF PTS
    cens.dist = c("EXP", "LN", "UNIF"),  ## CENSOR DIST.
    # meanlog, sdlog, minc, maxc, 
    marker.par
    # x1.u, x2.u, x1.s, x2.s, v.sd ## MARKER'S DIST
){

  # S <- 30; tK <- 2; n.mean <-4; n.sd <-1; Exp.r <- 0.2; cens.dist = "EXP"; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 0.2; N.min = 50; N.max = 150
  cens.dist <- match.arg(cens.dist)
  
  
  dt <- sapply(1:S, function(i) {
    
    ## 2. FOR i-th STUDY ----
    
    ## NUMBER OF SUBJECTS FOR S STUDIES: n0, n1
    
    n.i <- round(runif(1, min = npt.range[1], max = npt.range[2]))
    # n.i <- max(4, round(rlnorm(1, mean = n.mean, sd = n.sd)))
    # i <- 1; n.i <- 10
    
    ## 2.1 TIMES (YEARS) ----
    
    T.fail <- exp(rnorm(n.i, mean = 1, sd = 1))
    
    
    if (cens.dist == "EXP")  T.cens <- rexp(n.i, rate = 0.2)
    if (cens.dist == "LN")   T.cens <- rlnorm(n.i, meanlog = meanlog, sdlog = sdlog)
    if (cens.dist == "UNIF") T.cens <- runif(n.i, min = minc, max = maxc)
    
    Y.time <- pmin(T.fail, T.cens)
    status <- as.numeric(T.fail==Y.time)
    
    # MEDIAN YEARS
    
    med.ty <- median(Y.time)
    
    
    ## 2.2. BIOMARKERS ----
    
    x1.u <- marker.par[1] 
    x2.u <- marker.par[2] 
    x1.s <- marker.par[3] 
    x2.s <- marker.par[4] 
    v.sd <- marker.par[5] 
    X1 <- x1.u + x1.s*rlogis(n.i)
    X2 <- x2.u + x2.s*rlogis(n.i)
    X <- ifelse(T.fail <= tK, X1, X2)
    
    ## CUTOFF VALUES FOR i-th STUDY ----
    
    v.i <- rnorm(1, mean = 0.5, sd = v.sd)
    # v.i <- runif(1, range(X))
    ## 2.3. n0, n1 ----
    
    X.grp <- ifelse(X >= v.i, 1, 0)  #
    
    n1 <- sum(X.grp)          # HEG
    n0 <- n.i - n1  # LEG
    
    if(0 %in% c(n0, n1)) {
      
      x <- c(i, tK, n0, n1, rep(NA, 7), v.i)
      
      names(x) <- c("study.id", "ty", "n0", "n1", "s0", "s1", "med.ty", "s0_mct", "s1_mct", "u_lnHR", "v_lnHR", "cutoff")
      
      # fit.km <- NULL
      
      # list(s.tK = x, fit.km = NULL, marker = X)
      
    } else {
      
      fit.km <- survfit(Surv(Y.time, status) ~ X.grp)
      # km.info <- summary(fit.km)
      # plot(fit.km)
      
      ## 2.4. s1_s0 AT tK----
      t0 <- fit.km$time[1:n0]
      t1 <- fit.km$time[-c(1:n0)]
      
      surv0 <- fit.km$surv[1:n0]
      surv1 <- fit.km$surv[-c(1:n0)]
      
      if (sum(t0 <= tK) == 0) s0.by.tK <- 1 else s0.by.tK <- min(surv0[t0 <= tK])
      if (sum(t1 <= tK) == 0) s1.by.tK <- 1 else s1.by.tK <- min(surv1[t1 <= tK])
      
      ## 2.5. s1_mct, s0_mct AT med.ty----
      
      
      ## DIFFERENCE FROM Y.TIME
      
      if (sum(t0 <= med.ty) == 0) s0_mct <- 1 else s0_mct <- min(surv0[t0 <= med.ty])
      if (sum(t1 <= med.ty) == 0) s1_mct <- 1 else s1_mct <- min(surv1[t1 <= med.ty])
      
      
      ## 2.6. u_lnHR (FROM COX REGRESSION)
      if(testit::has_warning(coxph(Surv(Y.time, status) ~ X.grp))) u_lnHR <- v_lnHR <- NA else {
        
        fit.cox <- coxph(Surv(Y.time, status) ~ X.grp)
        
        u_lnHR  <- fit.cox$coefficients 
        v_lnHR  <- fit.cox$var
        se_lnHR <- sqrt(v_lnHR)
        # t_lnHR  <- u_lnHR / se_lnHR
        
      }
      
      x <- c(i, tK, n0, n1, s0.by.tK, s1.by.tK, med.ty, s0_mct, s1_mct, u_lnHR, v_lnHR, v.i)
      
      names(x) <- c("study.id", "ty", "n0", "n1", "s0", "s1", "med.ty", "s0_mct", "s1_mct", "u_lnHR", "v_lnHR", "cutoff")
      
      # plot(fit.km)
    }
    
    x
    # return(list(s.tK = x, fit.km = fit.km, marker = X))
    
  })
  
  dtt <- as.data.frame(t(dt))
  return(dtt)
  
}


.cens.eta <- function(data, med.year, n1, n0, s1_med, s0_med, eta.range = c(0,1), init.eta = 0.01){
  
  dt.mct <- data.frame(n1 = eval(substitute(n1), data), 
                       n0 = eval(substitute(n0), data),
                       med.year = eval(substitute(med.year), data),
                       s1_med = eval(substitute(s1_med), data),
                       s0_med = eval(substitute(s0_med), data))
  
  S_tf     <-  with(dt.mct, (s0_med * n0 + s1_med * n1) / (n0+n1) )
  
  U <- function(eta) sum( (0.5 - S_tf * exp(- dt.mct$med.year * eta))^2 )  ## G = P(C>t) =  exp(- eta t)
  
  int.eta <- c(eta = init.eta)
  
  nlminb(int.eta, U, lower = eta.range[1], upper = eta.range[2])
  
  
}

.convert.dt <- function(data, tK, study, ty, n1, n0, s1, s0, u_lnHR, v_lnHR, eta){
  
  # data <- p_dt.ost.df
  dt.os <- data.frame(study = eval(substitute(study), data),
                      ty =  eval(substitute(ty), data),
                      n1 = eval(substitute(n1), data), 
                      n0 = eval(substitute(n0), data),
                      s1 = eval(substitute(s1), data),
                      s0 = eval(substitute(s0), data),
                      u_lnHR = eval(substitute(u_lnHR), data), 
                      v_lnHR = eval(substitute(v_lnHR), data)
  )
  
  dt.os$t_lnHR <- with(dt.os, u_lnHR / sqrt(v_lnHR))
  uniq.study <- sort(unique(dt.os$study))
  dt.ost <- dt.os[dt.os$ty == tK,  ] 
  ## CALCULATE: S1, S0, q1, q0
  
  S1 <- dt.ost$s1
  S0 <- dt.ost$s0
  q1  <- with(dt.ost, n1 / (n0+n1))											
  q0  <- with(dt.ost, n0 / (n0+n1))				
  
  
  ## AD_HOC CORRECTION
  
  S1 <- ifelse((S1 == 1), with(dt.ost, (n1 + 0.5) / (n1 + 1)), S1)
  S1 <- ifelse((S1 == 0), with(dt.ost,       0.5  / (n1 + 1)), S1)
  S0 <- ifelse((S0 == 1), with(dt.ost, (n0 + 0.5) / (n0 + 1)), S0)
  S0 <- ifelse((S0 == 0), with(dt.ost,       0.5  / (n0 + 1)), S0)
  
  
  # x <- S1
  # y <- S0
  # z <- q1
  # w <- q0
  
  
  ## logit-SE & logit-SP
  
  sen <- (1-S1) * q1 / ( (1-S1) * q1 + (1-S0) * q0 )
  spe <- S0 * q0 / (S1 * q1 + S0 * q0)
  
  dt.ost$u_sen <- qlogis(sen)
  dt.ost$u_spe <- qlogis(spe)
  
  ## GREENWOOD VARIANCE
  
  green.var <- vapply(uniq.study, function(i){
    
    # i <- 1
    study.i   <- dt.os[dt.os$study == i,]
    study.i$G <- exp(- study.i$ty * eta)
    study.i   <- rbind.data.frame(rep(1, ncol(study.i)), study.i)
    
    pos.t        <- match(tK, study.i$ty)
    
    if(!is.na(pos.t)){
      
      ## sigma0^2
      
      study.i$inv.S0G  <- 1 / (study.i$s0^2 * study.i$G)
      study.i$half.h0  <- c(0, diff(study.i$s0, lag = 1))/2
      study.i$lag.sum0 <- c(0, head(study.i$inv.S0G, -1) + tail(study.i$inv.S0G, -1))
      study.i$trapez0  <- study.i$lag.sum0 * study.i$half.h0
      study.i$cum.sum0 <- cumsum(study.i$trapez0)
      
      sigma0_2 <- -study.i[study.i$ty == tK, "s0"]^2 * study.i$cum.sum0[pos.t]
      
      
      ## sigma1^2
      
      study.i$inv.S1G  <- 1/(study.i$s1^2 *study.i$G)
      study.i$half.h1  <- c(0, diff(study.i$s1, lag = 1))/2
      study.i$lag.sum1 <- c(0,head(study.i$inv.S1G, -1) + tail(study.i$inv.S1G, -1))
      study.i$trapez1  <- study.i$lag.sum1 * study.i$half.h1
      study.i$cum.sum1 <- cumsum(study.i$trapez1)
      
      sigma1_2 <- -study.i[study.i$ty == tK, "s1"]^2 * study.i$cum.sum1[pos.t]
      
      
      c(i, sigma1_2, sigma0_2)
      
    } else {c(i, NA, NA)}
    
  }, c("study" = 0,"km.v1" = 0, "km.v0" = 0))
  
  km.var.mat  <- t(green.var)
  km.var.comp <- km.var.mat[complete.cases(km.var.mat), ]
  
  dt.ost <- merge(dt.ost, km.var.comp, by = "study" , all.x = TRUE)
  
  
  ##
  ## 2.3. H/n (2x2) MATRIX ----
  ##
  
  ## PARTIAL DERIVATIVES
  
  g_senx <-  1 / (S1-1)
  g_seny <-  1 / (1-S0)
  g_senz <-  1 / q1
  g_senw <- -1 / q0
  
  g_spex <- -1 / S1
  g_spey <-  1 / S0
  g_spez <- -1 / q1
  g_spew <-  1 / q0
  
  inv.q1 <-  1 / q1
  inv.q0 <-  1 / q0
  q1.q0  <-  q1*q0
  
  # km.v1 <- dt.ost$km.v1
  # v0 <- dt.ost$km.v0
  n  <- dt.ost$n1 + dt.ost$n0
  
  dt.ost$v_sen <- with(dt.ost, g_senx^2 * inv.q1 * km.v1 + g_seny^2 * inv.q0 * km.v0 + (g_senz - g_senw)^2 * q1.q0) / n
  dt.ost$v_spe <- with(dt.ost, g_spex^2 * inv.q1 * km.v1 + g_spey^2 * inv.q0 * km.v0 + (g_spez - g_spew)^2 * q1.q0) / n
  dt.ost$v_senspe  <- with(dt.ost, g_senx * g_spex * inv.q1 * km.v1 + g_seny * g_spey * inv.q0 * km.v0 + (g_senz - g_senw) * (g_spez - g_spew) * q1.q0) / n
  
  ##
  ## 2.4. FOR STUDY I, CALCULATE: cov(R1, W), cov(R0, W)   ----
  ##
  
  inv.B <- dt.ost$v_lnHR
  
  ## cov(R, W) FOR ALL STUDIES
  
  int <- vapply(uniq.study, function(i) {
    
    # i <- 1
    
    study.i <- dt.os[dt.os$study == i,]
    n1 <- study.i$n1[1]
    n0 <- study.i$n0[1]
    
    study.i <- rbind(rep(1, ncol(study.i)), study.i)
    
    pos.t <- match(tK, study.i$ty)
    
    if(!is.na(pos.t)){
      
      study.i$inv.S    <- with(study.i, (n0 + n1) / (n0 * s0 + n1 * s1))
      study.i$lag.sum  <- c(0, head(study.i$inv.S, -1) + tail(study.i$inv.S, -1))
      study.i$half.h0  <- c(0, diff(study.i$s0, lag = 1))/2
      study.i$half.h1  <- c(0, diff(study.i$s1, lag = 1))/2
      study.i$trapez0  <- study.i$lag.sum * study.i$half.h0
      study.i$trapez1  <- study.i$lag.sum * study.i$half.h1
      study.i$cum.sum0 <- cumsum(study.i$trapez0)
      study.i$cum.sum1 <- cumsum(study.i$trapez1)
      int0             <- study.i[pos.t, "cum.sum0"]
      int1             <- study.i[pos.t, "cum.sum1"]
      
      c(i, int1, int0)
      
    } else c(i, NA, NA)
    
  }, c("study" = 0,"int1" = 0, "int0" = 0))
  
  int.mat <- t(int)
  int.mat.comp <- int.mat[complete.cases(int.mat), ]
  
  dt.ost <- merge(dt.ost, int.mat.comp, by = "study" , all.x = TRUE)
  
  dt.ost$v_senlnHR <- with( dt.ost, g_senx * S1 * inv.q1 * inv.B * (log(S1) - int1) + g_seny * S0 * inv.q0 * inv.B * (log(S0) - int0) ) /n 
  dt.ost$v_spelnHR <- with( dt.ost, g_spex * S1 * inv.q1 * inv.B * (log(S1) - int1) + g_spey * S0 * inv.q0 * inv.B * (log(S0) - int0) ) /n
  
  dt.ost <-  subset(dt.ost, select = -c(km.v1, km.v0, int1, int0, n1, n0, s1, s0) )
  
  names(dt.ost) <- c("study", "tk","y3","v3","t","y1","y2","v1","v2","v12","v13","v23")
  
  return(dt.ost[c(1,2,6,7,3,8,9,4,10:12)])
  
  
}

##
## GENERATE LOGIT PDATA AND SDATA
##
.ps.sim.logit <- function(
    S, 
    tK, 
    npt.range,  ## NUMBER OF PTS
    cens.dist = c("EXP", "LN", "UNIF"),  ## CENSOR DIST.
    # meanlog, sdlog, minc, maxc, 
    marker.par,
    alpha.n,
    beta
    # x1.u, x2.u, x1.s, x2.s, v.sd ## MARKER'S DIST
){
  # S=20; npt.range=set.ls$npt.range1; marker.par=set.ls$marker.par1; alpha.n=1; tK=2
  p_dataSt <- .p.sim.ipd(
    S, 
    tK, 
    npt.range,  ## NUMBER OF PTS
    cens.dist = cens.dist,  ## CENSOR DIST.
    # meanlog, sdlog, minc, maxc, 
    marker.par)%>%na.omit()
  
  eta <- .cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
  dataSt.p <- .convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)%>%na.omit()
  
  t = dataSt.p$y3/sqrt(dataSt.p$v3)
  p.select <- pnorm(alpha[alpha.n] + beta*t)
  # round(p.select,3)
  select <- rbinom(nrow(dataSt.p),1, p.select)
  # c(mean(p.select[select==1]),mean(p.select[select==0]))
  s_dataSt <- p_dataSt[select==1,]
  eta.s <- .cens.eta(data = s_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
  dataSt.s <- .convert.dt(data = s_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta.s)%>%na.omit()
  
  data = list(pdata = dataSt.p, sdata = dataSt.s)
  return(data)
}

##
## GENERATE LOGIT PDATA ONLY
##
.p.sim.logit <- function(
    S, 
    tK, 
    npt.range,  ## NUMBER OF PTS
    cens.dist = c("EXP", "LN", "UNIF"),  ## CENSOR DIST.
    # meanlog, sdlog, minc, maxc, 
    marker.par
    # x1.u, x2.u, x1.s, x2.s, v.sd ## MARKER'S DIST
){
  # S=20; npt.range=set.ls$npt.range1; marker.par=set.ls$marker.par1; alpha.n=1; tK=2
  p_dataSt <- .p.sim.ipd(
    S, 
    tK, 
    npt.range,  ## NUMBER OF PTS
    cens.dist = cens.dist,  ## CENSOR DIST.
    # meanlog, sdlog, minc, maxc, 
    marker.par)%>%na.omit()
  
  eta <- .cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
  dataSt.p <- .convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)%>%na.omit()
  
  pdata = dataSt.p
  return(pdata)
}







