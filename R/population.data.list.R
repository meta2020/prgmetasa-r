##
## 1. CREATE POPULATION STUDIES DATA (LIST) ----
##

population.data.list <- function(
	S, 
	tK, 
	N.min, N.max,  ## NUMBER OF PTS
	CD = c("EXP", "LN", "UNIF"),  ## CENSOR DIST.
	Exp.r, meanlog, sdlog, minc, maxc, 
	x1.u, x2.u, x1.s, x2.s, v.sd
	){
	
	# S <- 30; tK <- 2; n.mean <-4; n.sd <-1; Exp.r <- 0.2; CD = "EXP"; x1.u <- 0.7; x2.u <- 0.3; x1.s <- 0.1; x2.s <- 0.3; v.sd <- 1; N.min = 50; N.max = 150
	
	lapply(1:S, function(i){
		
		## 2. FOR i-th STUDY ----
		
		## NUMBER OF SUBJECTS FOR S STUDIES: n0, n1
		
		n.i <- round(runif(1, min = N.min, max = N.max))
	  # n.i <- max(4, round(rlnorm(1, mean = n.mean, sd = n.sd)))
	  # i <- 1; n.i <- 10
		
		## 2.1 TIMES (YEARS) ----
		
		T.fail <- exp(rnorm(n.i, mean = 1, sd = 1))
		
		if (CD == "EXP")  T.cens <- rexp(  n.i, rate = Exp.r)
		if (CD == "LN")   T.cens <- rlnorm(n.i, meanlog = meanlog, sdlog = sdlog)
		if (CD == "UNIF") T.cens <- runif( n.i, min = minc, max = maxc)
		
		Y.time <- pmin(T.fail, T.cens)
		status <- as.numeric(T.fail==Y.time)
		
		# MEDIAN YEARS
		
		med.ty <- median(Y.time)
		
		
		## 2.2. BIOMARKERS ----
		
		
		X1 <- x1.u + x1.s*rlogis(n.i)
		X2 <- x2.u + x2.s*rlogis(n.i)
		X <- ifelse(T.fail <= 2, X1, X2)
		
		## CUTOFF VALUES FOR i-th STUDY ----
		
		v.i <- rnorm(1, mean = 0.5, sd = v.sd)
		
		## 2.3. n0, n1 ----
		
		X.grp <- ifelse(X >= v.i, 1, 0)  #
		
		n1 <- sum(X.grp)          # HEG
		n0 <- n.i - n1  # LEG
		
		if(0 %in% c(n0, n1)) {
			
			x <- c(i, tK, n0, n1, rep(NA, 7), v.i)
			
			names(x) <- c("study.id", "ty", "n0", "n1", "s0", "s1", "med.ty", "s0_mct", "s1_mct", "u_lnHR", "v_lnHR", "cutoff")
			
			fit.km <- NULL
			
			# list(s.tK = x, fit.km = NULL, marker = X)
			
		} else{
			
			fit.km <- survfit(Surv(Y.time, status) ~ X.grp)
			# km.info <- summary(fit.km)
			plot(fit.km)
			
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
		
			return(list(s.tK = x, fit.km = fit.km, marker = X))
			
		
		
	})
	
}