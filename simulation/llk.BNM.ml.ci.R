##******************************************************************************
##
## LIKELIHOOD OF THE OBSERVED (NEGATIVE)
##
##******************************************************************************

.llk.BNM.ml <- function(par,y1, y2, v1, v2, v12){
  
  u1 <- par[1]
  u2 <- par[2]
  
  t1 <- par[3]
  t2 <- par[4]
  
  r  <- par[5]
  
  t11 <- t1^2
  t22 <- t2^2
  t12 <- t1*t2*r
  
  v11 <- v1  + t11
  v22 <- v2  + t22
  v12 <- v12 + t12
  
  
  ##
  ##  LOGLIKELIHOOD-1 OF y|Sigma ----
  ##
  
  det.vec <- v11 * v22 -v12^2
  
  log.det.vec <- suppressWarnings(log(det.vec))
  
  f.l1  <- ((y1-u1)^2*(v22) - 2*(y2-u2)*(y1-u1)*v12 + (y2-u2)^2*(v11)) / det.vec + log.det.vec
  
  s.l1  <- 0.5*sum(f.l1, na.rm = TRUE)
  
  return(s.l1)
}

BNM.mle <- function(y1, y2, v1, v2, v12){
  
  eps <- 0.001
  fn  <- function(par) .llk.BNM.ml(par, y1, y2, v1, v2, v12)

  fn.est <- try(nlminb(
    runif(5,-0.1,0.1),  ## initial values for mu1,2 tau1,2 rho1
    fn,
    lower = c(-5,-5,eps,eps,-0.999),
    upper = c(5,5,2,2,-eps)
  ))

  
  if(!inherits(fn.est,"try-error")) {
    
    num.hessian <- hessian(fn, fn.est$par)
    rownames(num.hessian) <- colnames(num.hessian) <- c("mu1","mu2","tau1","tau2","rho1")
    
    fn.est$par <- c(fn.est$par[1:2], NA, fn.est$par[3:4], NA, fn.est$par[5], NA, NA)
    names(fn.est$par) <- c("mu1","mu2","mu3","tau1","tau2","tau3","rho1","rho2","rho3")
    
    
    var.ml <- try(solve(num.hessian))
    
    if(!inherits(var.ml,"try-error")) {
      
      var.matrix <-  var.ml
      sauc.ci <- .SAUC.ci(fn.est$par[c(1,2,4,5,7)], var.matrix = var.matrix)
      
    } else {
      
      sauc <- .SAUC(fn.est$par[c(1,2,4,5,7)])
      sauc.ci = c(sauc=sauc, sauc.lb=NA, sauc.ub=NA)
      
    }
    
    res  <- c(fn.est, sauc.ci)
    
  } else res  <- NULL
  
  
  return(res)
  
}

# BNM.mle(pdata$y1, pdata$y2, pdata$v1, pdata$v2,pdata$v12)
