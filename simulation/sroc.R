.SROC <- function(
    par,  ## mu1 mu2 tau1 tau2 rho
    add = FALSE,
    ncols = NULL,
    add.spoint=TRUE,
    spoint.pch = 18,
    spoint.cex = 2,
    xlab = "FPR",
    ylab = "TPR",
    ...
){
  
  par <- as.matrix(par)
  
  if (!add) plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab = xlab, ylab = ylab, ...)
  
  if (is.null(ncols)) ncols <- gray.colors(ncol(par), gamma = 1, start = 0, end = 0.8)
  
  for (i in 1:ncol(par)) {
    
    u1  <- par[1,i]
    u2  <- par[2,i]
    t1  <- par[3,i]
    t2  <- par[4,i]
    r   <- par[5,i]
    
    roc <- function(x) plogis(u1 - (t1*r/t2) * (qlogis(x) + u2))
    curve(roc, 0, 1, col = ncols[i], add = TRUE,
          lty = 1, lwd = 1)
  }
  
  
  if (add.spoint) {
    sens <- plogis(par[1,])
    spec <- plogis(par[2,])
    points(1-spec, sens, col=ncols, pch = spoint.pch, cex = spoint.cex)
  }
  
}


##
## PLOT THE SROC CURVE USING LOGIT DATA
.SROC.reitsma <- function(data){
  
  fit =  BNM.mle(data$y1, data$y2, data$v1, data$v2, data$v12)
  par = fit$par
  
  plot(plogis(-data$y2),plogis(data$y1), xlim=c(0,1), ylim=c(0,1), xlab = "1-sp", ylab = "se", pch=4, col="grey")
  abline(a=0, b=1, col="grey")
  u1  <- par[1]
  u2  <- par[2]
  t1  <- par[4]
  t2  <- par[5]
  r   <- par[7]
  
  roc <- function(x) plogis(u1 - (t1*r/t2) * (qlogis(x) + u2))
  curve(roc, 0, 1, add = TRUE, lty = 1, lwd = 1)
  legend("bottomright", legend=round(c(fit$sauc*100, r),5))
  
  sens <- plogis(par[1])
  spec <- plogis(par[2])
  points(1-spec, sens, pch = 18, cex=2)
  

}
.SROC.reitsma.add <- function(data){
  
  fit =  BNM.mle(data$y1, data$y2, data$v1, data$v2, data$v12)
  par = fit$par
  points(plogis(-data$y2),plogis(data$y1), pch=4, col="red")
  
  u1  <- par[1]
  u2  <- par[2]
  t1  <- par[4]
  t2  <- par[5]
  r   <- par[7]
  
  roc <- function(x) plogis(u1 - (t1*r/t2) * (qlogis(x) + u2))
  curve(roc, 0, 1, add = TRUE, col=2)
  legend("topright", legend=round(c(fit$sauc*100, r),5), col="red")
  
  sens <- plogis(par[1])
  spec <- plogis(par[2])
  points(1-spec, sens, pch = 18, cex=2, col=2)
  
  
}


.DID.sroc <- function(u1, u2, t1, t2, r, var.matrix){
  
  Q1 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    g*(1-g)
  }
  
  Q2 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.u2 <- (-r*t1/t2)*g*(1-g)
  }
  
  Q3 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.t1 <- (-r/t2*(qlogis(x)+u2))*g*(1-g)
    
  }
  
  Q4 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.t2 <- r*t1/t2^2*( qlogis(x)+u2)*g*(1-g)
    
  }
  
  Q5 <- function(x) {
    
    g <- plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
    p.r  <- (-t1)/t2*(qlogis(x) + u2)*g*(1-g)
  }
  
  fd <- c(integrate(Q1, 0, 1)$value,
          integrate(Q2, 0, 1)$value,
          integrate(Q3, 0, 1)$value,
          integrate(Q4, 0, 1)$value,
          integrate(Q5, 0, 1)$value
  )
  
  (fd %*% var.matrix %*% fd)
  
}


.SAUC.ci <- function(
    par, 
    var.matrix, 
    ci.level = 0.95){
  
  u1 <- par[1]; u2 <- par[2]
  t1 <- par[3]; t2 <- par[4]; r <- par[5]
  

  sroc <- function(x) plogis(u1 - (t1*t2*r/(t2^2)) * (qlogis(x) + u2))
  sauc.try <- try(integrate(sroc, 0, 1))
  if(!inherits(sauc.try, "try-error")) sauc <- sauc.try$value else sauc <- NA
  
  sauc.lb <-  plogis(qlogis(sauc) + qnorm((1-ci.level)/2, lower.tail = TRUE) *
                       suppressWarnings(
                         sqrt(.DID.sroc(u1, u2, t1, t2, r, var.matrix))/(sauc*(1-sauc))) )
  
  sauc.ub <-  plogis(qlogis(sauc) + qnorm((1-ci.level)/2, lower.tail = FALSE)*
                       suppressWarnings(
                         sqrt(.DID.sroc(u1, u2, t1, t2, r, var.matrix))/(sauc*(1-sauc))) )
    
  
  sauc.ci <- c(sauc, sauc.lb, sauc.ub)
  names(sauc.ci) <- c("sauc", "sauc.lb", "sauc.ub")
  
  return(sauc.ci)
}

