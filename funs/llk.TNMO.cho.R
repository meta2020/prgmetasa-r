##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.TNMO.cho <- function(par,
                    y1, y2, y3,
                    s11, s22, s33,
                    s12, s13, s23,
                    p,
                    a.interval
                  ){

  
  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  t1  <- par[4]     
  t2  <- par[5]    
  t3  <- par[6]
  
  t4  <- par[7]      
  t5  <- par[8]
  t6  <- par[9]      
  
  b   <- par[10]
  
  L <- diag(rep(0,3))
  L[lower.tri(L, diag = TRUE)] <- c(t1, t2, t3, t4, t5, t6)
  Omega <- L %*% t(L)
  
  Sigma <- cbind(s11, s12, s13, s12, s22, s23, s13, s23, s33)
  
  t_lnHR   <- y3/sqrt(s33)

  ##
  ## FUNCTOIN b(Sigma) ----
  ##

  t33 <- (t3^2 + t5^2 + t6^2)

  f.b <- function(a){

    sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / s33)))

    pnorm( (a + b * u3 / sqrt(s33)) / sq )

  }


  ##
  ## FIND THE ROOT OF a = a.opt ----
  ##

  a.p <- function(a) mean(1/f.b(a), na.rm = TRUE) - 1/p

  a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="yes"), silent = FALSE))

  a.opt <- a.opt.try$root


  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  det <- sapply(1:nrow(Sigma), function(i) {
    
    det(matrix(Sigma[i, ], nrow = 3, ncol = 3) + Omega)
    
  }) 
  
  
  log.dev <- suppressWarnings(log(det)) 
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  
  Y <- cbind(Y1, Y2, Y3)
  
  uOu <- sapply(1:nrow(Sigma), function(i) {
    
    inv.try <- try(solve(matrix(Sigma[i, ], 3, 3) + Omega), silent = TRUE)
    if(!inherits(inv.try, "try-error")) inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega) else inv.det <- matrix(NA, 3, 3)
    
    # inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega)
    Y[i, ] %*% inv.det %*% Y[i, ]
    
  })  
  
  s.l1 <- -0.5 * (sum(log.dev + uOu, na.rm = TRUE))
  
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
 # -s.l1
  
}
