##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************


llk.BNM.cho <- function(par, y1, y2, v1, v2, v12){

  u1  <- par[1]      ## logit-se
  u2  <- par[2]      ## logit-sp
  
  t1  <- par[3]      ## logit-se
  t2  <- par[4]    
  t3  <- par[5]
  
  
  L <- diag(rep(0,2))
  L[lower.tri(L, diag = TRUE)] <- c(t1, t2, t3)
  Omega <- L %*% t(L)
  
  Sigma <- cbind(v1, v12, v12, v2)
    
  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##

  det <- sapply(1:nrow(Sigma), function(i) {
    
    det(matrix(Sigma[i, ], nrow = 2, ncol = 2) + Omega)
    
  }) 
  
  
  log.dev <- suppressWarnings(log(det)) 
  
  Y1  <- y1 - u1
  Y2  <- y2 - u2
  
  Y   <- cbind(Y1, Y2) 
  
  uOu <- sapply(1:nrow(Sigma), function(i) {
    
    inv.try <- try(solve(matrix(Sigma[i, ], 2, 2) + Omega), silent = TRUE)
    if(!inherits(inv.try, "try-error")) inv.det <- solve(matrix(Sigma[i, ], 2, 2) + Omega) else inv.det <- matrix(NA, 2, 2)

    # inv.det <- solve(matrix(Sigma[i, ], 3, 3) + Omega)
    Y[i, ] %*% inv.det %*% Y[i, ]
    
  })  
  
  s.l1 <- -0.5 * (sum(log.dev + uOu, na.rm = TRUE))
  
  -s.l1
}
