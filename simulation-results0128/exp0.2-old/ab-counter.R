

## POPULATION DATA
S <- 50
tK <- 2
N.min <- 50
N.max <- 150
CD <- "EXP"
Exp.r <- 0.2
x1.u <- 0.7
x2.u <- 0.3 
x1.s <- 0.1
x2.s <- 0.3
v.sd <- 0.2

p_dataSt.all <- population.data.list(S, tK, N.min, N.max, CD, Exp.r, meanlog, sdlog, minc, maxc, x1.u, x2.u, x1.s, x2.s, v.sd)
p_dataSt <- as.data.frame(t(sapply(1:length(p_dataSt.all), function(i) p_dataSt.all[[i]][[1]])))

## REMOVE NA

p_dataSt <- p_dataSt[complete.cases(p_dataSt),]

## ANALYTICAL POPULATION DATA 
eta <- cens.eta(data = p_dataSt, med.year = med.ty, n1 = n1, n0 = n0, s1_med = s1_mct, s0_med = s0_mct)$par
dataSt.p <- convert.dt(data = p_dataSt, tK, study = study.id, ty = ty, n1 = n1, n0 = n0, s1 = s1, s0 = s0, u_lnHR = u_lnHR, v_lnHR = v_lnHR, eta = eta)

## REMOVE NA

dataSt.p <- dataSt.p[complete.cases(dataSt.p),]

## 2. CREATE OBSERVED LIST ----




b <- seq(0,5,0.5)
atb <- seq(min(b),max(b),0.5)
a <- seq(-7,7,0.5)
ata <- seq(min(a),max(a),0.5)

ba <- expand.grid(b=b,a=a)
n <- nrow(ba)
lb <- length(b)
la <- length(a)


bap.list <- lapply(1:1, function(i){
  
  p.a.e <- vapply(1:n, function(j) {
    
    t    <- dataSt.p$t_lnHR
    mean <- mean(pnorm(ba$a[j] + ba$b[j]*t, lower.tail = TRUE))
    med  <- median(pnorm(ba$a[j] + ba$b[j]*t, lower.tail = TRUE))
    
    c(mean, med)
    
  },  c(mean = 0, median = 0))
  
  colnames(p.a.e) <- paste0(ba$b, ";", ba$a)
  
  bap <- matrix(p.a.e[1,], lb, la, dimnames = list(b, a))
  bap2 <- matrix(p.a.e[2,], lb, la, dimnames = list(b, a))
  
  list(bap.e = bap,
       bap.m = bap2,
       p.a.e = p.a.e)
  
})

ab.contour <- function(i){
  
  z  <- t(bap.list[[i]]$bap.e)
  z2 <- t(bap.list[[i]]$bap.m)
  pe <- bap.list[[i]]$p.a.e
  #pos<- which(round(ba$b,2) == by & round(ba$a,2) == ax)
  
  # contour(a, b, z,  nlevels = 10, axes=FALSE, frame.plot = FALSE)
  contour(a, b, z2,  nlevels = 10, axes=FALSE, add = TRUE, col = 2, lty = 3)
  #points(ax, by, col=2, pch=4)
  axis(1, at = ata)
  axis(2, at = atb)
  grid()
  # title(main = paste(sprintf("%s = %.1f", c("t1","t2","r"), c(X[3:5,i])), collapse="; "), xlab="a", ylab = "b")
  
}

 
ab.contour(1)

forest(m1, 
       layout = "RevMan5", common = FALSE, random = FALSE,
       hetstat = FALSE,
       
               label.right = "Favours control",
               label.left = "Favours experimental",
               prediction = FALSE)
