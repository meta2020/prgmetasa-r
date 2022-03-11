##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.TNMO.mat <- function(par,
                    y1, y2, y3,
                    s11, s22, s33,
                    s12, s13, s23,
                    p,
                    a.interval
                  ){

  
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
    
    sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / s33)))
    
    pnorm( (a + b * u3/sqrt(s33)) / sq )
    
  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##

  a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}

  a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="yes"), silent = TRUE)) 

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  y.hat <- cbind(y1, y2, y3)
  V_O   <- cbind(v11, v12, v13, v12, v22, v23, v13, v23, v33)
  
  det <- sapply(1:nrow(V_O), function(i) {
    
    det(matrix(V_O[i, ], 3, 3))
    
  }) 
  
  
  log.dev <- suppressWarnings(log(det)) 
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  
  Y <- cbind(Y1, Y2, Y3)
  
  uOu <- sapply(1:nrow(V_O), function(i) {
    
    inv.det <- solve(matrix(V_O[i, ], 3, 3))
    Y[i, ] %*% inv.det %*% Y[i, ]
    
  })  
  
  f.l1 <- log.dev + uOu
  
  s.l1 <- -0.5 * sum(f.l1, na.rm = TRUE)


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
