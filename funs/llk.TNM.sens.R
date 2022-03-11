##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.TNM.o <- function(par,
                    y1, y2, y3,
                    s11, s22, s33,
                    s12, s13, s23,
                    p,
                    a.interval
                  ){
# par <- par0

  n   <- length(y1) 
  
  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  t1  <- par[4]      ## logit-se
  t2  <- par[5]    
  t3  <- par[6]
  
  r1  <- par[7]      ## logit-sp
  r2  <- par[8]
  r3  <- par[9]      ## ln-HR
  
  b   <- par[10] 

  t11 <- t1^2
  t22 <- t2^2
  t33 <- t3^2

  t12 <- t1*t2*r1
  t23 <- t2*t3*r2
  t13 <- t1*t3*r3
  
  # y1 = dataSt$u_sen
  # y2 = dataSt$u_spe
  # y3 = dataSt$u_lnHR
  # 
  # v1 = dataSt$v_sen
  # v2 = dataSt$v_spe
  # v3 = dataSt$v_lnHR
  # 
  # v12 = dataSt$v_senspe
  # v13 = dataSt$v_senlnHR
  # v23 = dataSt$v_spelnHR
  # 
  # u1 <- model3$coefficients[1]
  # u2 <- model3$coefficients[2]
  # u3 <- model3$coefficients[3]
  # 
  # t11 <- model3$Psi[1,1]
  # t22 <- model3$Psi[2,2]
  # t33 <- model3$Psi[3,3]
  # 
  # t12 <- model3$Psi[1,2]
  # t13 <- model3$Psi[1,3]
  # t23 <- model3$Psi[2,3]
  
  
  v11 <- s11  + t11
  v22 <- s22  + t22
  v33 <- s33  + t33
  
  v12 <- s12 + t12 
  v13 <- s13 + t13
  v23 <- s23 + t23
  

  t_lnHR   <- y3/sqrt(s33)

  ##
  ## FUNCTOIN b(Sigma) ----
  ##

  f.b <- function(a){
    
    sq <- suppressWarnings(sqrt(1 + b^2 * (v33 / s33)))
    
    pnorm( (a + b * u3/sqrt(s33)) / sq )
    
  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##

  a.p <- function(a) {sum(1/f.b(a), na.rm = TRUE) - n/p}

  a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="yes"), silent = TRUE)) 

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  
  # sapply(1:nrow(V_O), function(i) {
  # 
  #   det(matrix(V_O[i, ], 3, 3))
  # 
  # }) - det
  
  det1  <- v22 * v33 - v23 * v23
  det2  <- v12 * v33 - v13 * v23
  det3  <- v12 * v23 - v13 * v22

  det   <- v11 * det1 - v12 * det2 + v13 * det3
  
  log.dev <- suppressWarnings(log(det)) 
  
  a11 <-  (v22 * v33 - v23 * v23)
  a12 <- -(v12 * v33 - v13 * v23)
  a13 <-  (v12 * v23 - v13 * v22)
  
  a21 <- -(v12 * v33 - v23 * v13)
  a22 <-  (v11 * v33 - v13 * v13)
  a23 <- -(v11 * v23 - v12 * v13)
  
  a31 <-  (v12 * v23 - v13 * v22)
  a32 <- -(v11 * v23 - v12 * v13)
  a33 <-  (v11 * v22 - v12 * v12)
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  
  y_Adj_y <- Y1 * (Y1 * a11 + Y2 * a21 + Y3 * a31) + 
    Y2 * (Y1 * a12 + Y2 * a22 + Y3 * a32) +
    Y3 * (Y1 * a13 + Y2 * a23 + Y3 * a33)
  
  # Y <- cbind(Y1, Y2, Y3)
  
  # sapply(1:nrow(V_O), function(i) {
  #   
  #   inv.det <- solve(matrix(V_O[i, ], 3, 3))
  #   Y[i, ] %*% inv.det %*% Y[i, ]
  #   
  # })  - y_inv.det_y / det
  
  f.l1 <- log.dev + y_Adj_y / det
  
  s.l1 <- -0.5 * sum(f.l1, na.rm = TRUE)
    
  # model3$logLik - (-0.5 * sum(f.l1, na.rm = TRUE))
  
  
  



  ##
  ##  LOGLIKELIHOOD-2 OF a(a.opt) ----
  ##


  f.l2 <- pnorm(a.opt + b * t_lnHR)

  s.l2 <- sum( log(f.l2), na.rm = TRUE )


  ##
  ##  LOGLIKELIHOOD-3 OF b(a.opt) ----
  ##

  f.l3 <- f.b(a.opt)

  s.l3 <- sum( log(f.l3), na.rm = TRUE )


  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  -(s.l1 + s.l2 - s.l3) ## NEGATIVE

  
}
