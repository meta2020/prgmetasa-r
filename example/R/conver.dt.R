convert.dt <- function(data, tK, study, ty, n1, n0, s1, s0, u_lnHR, v_lnHR, eta){
	
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
	
	dt.ost <-  subset(dt.ost, select = -c(km.v1, km.v0, int1, int0) )
	
	dt.ost
	
}