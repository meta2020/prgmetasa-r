##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


.llk.TNM.cho <- function(par, y1, y2, y3, v1, v2, v3, v12, v13, v23 ){

  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  u3  <- par[3]      ## ln-HR
  
  t1  <- par[4]      ## logit-se
  t2  <- par[5]    
  t3  <- par[6]
  
  t4  <- par[7]      ## logit-sp
  t5  <- par[8]
  t6  <- par[9]      ## ln-HR
  
  L <- diag(rep(0,3))
  L[lower.tri(L, diag = TRUE)] <- c(t1, t2, t3, t4, t5, t6)
  Omega <- L %*% t(L)
  
  Sigma <- cbind(v1, v12, v13, v12, v2, v23, v13, v23, v3)

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
  # 
  Y <- cbind(Y1, Y2, Y3) 
  
  uOu <- sapply(1:nrow(Sigma), function(i) {
    
    inv.try <- try(solve(matrix(Sigma[i, ], 3, 3) + Omega), silent = TRUE)
    if(!inherits(inv.try, "try-error")) inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega) else inv.det <- matrix(NA, 3, 3)

    Y[i, ] %*% inv.det %*% Y[i, ]
    
  })  
  
  s.l1 <- -0.5 * (sum(log.dev + uOu, na.rm = TRUE))
  
  ##
  ##  FINAL LOGLIKELIHOOD ----
  ##

  -s.l1
}
