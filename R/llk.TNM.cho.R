##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.TNM.cho <- function(par,
                    y1, y2, y3,
                    s11, s22, s33,
                    s12, s13, s23
                  ){

  # u1 <- u2 <- u3 <- 0.5; t1 <- t2 <- t3 <- 1; t4 <- t5 <- t6 <- 2; b <- 1
  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  # u <- c(u1, u2, u3)
  
  t1  <- par[4]     
  t2  <- par[5]    
  t3  <- par[6]
  
  t4  <- par[7]      
  t5  <- par[8]
  t6  <- par[9]      
  
  L <- diag(rep(0,3))
  L[lower.tri(L, diag = TRUE)] <- c(t1, t2, t3, t4, t5, t6)
  Omega <- L %*% t(L)
  
  Sigma <- cbind(s11, s12, s13, s12, s22, s23, s13, s23, s33)
    



  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  # y.hat <- cbind(y1, y2, y3)

  det <- sapply(1:nrow(Sigma), function(i) {
    
    det(matrix(Sigma[i, ], nrow = 3, ncol = 3) + Omega)
    
  }) 
  
  
  log.dev <- suppressWarnings(log(det)) 
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  Y3  <- y3 - u3
  # 
  Y <- cbind(Y1, Y2, Y3) 
  
  uOu <- sapply(1:nrow(Sigma), function(i) {
    
    inv.try <- try(solve(matrix(Sigma[i, ], 3, 3) + Omega), silent = TRUE)
    if(!inherits(inv.try, "try-error")) inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega) else inv.det <- matrix(NA, 3, 3)

    # inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega)
    Y[i, ] %*% inv.det %*% Y[i, ]
    
  })  
  
  s.l1 <- -0.5 * (sum(log.dev + uOu, na.rm = TRUE))
  

  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  # -(s.l1 + s.l2 - s.l3) ## NEGATIVE

  -s.l1
}
