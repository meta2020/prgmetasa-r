##
## FUNCTION USED IN THE EXAMPLE
##


progmeta.sa <- function(
    y1, y2, y3, v1, v2, v3, v12, v13, v23, 
    p, 
    sauc.type = "sroc",
    a.interval= c(-10, 10)
    ){

p <- p

fn.o <- function(par) clk.TNM.ml(par, y1, y2, y3, v1, v2, v3, v12, v13, v23,
                                  p = p, 
                                  a.interval= a.interval
)

eps <- sqrt(.Machine$double.eps)



opt.o <- nlminb(c(rep(0,3), rep(0.1, 3), rep(-0.1, 3), 0.5), fn.o,
                lower = c(rep(-Inf,3), rep(0,3), rep(-1,3), 0 ),
                upper = c(rep( Inf,3), rep(Inf,3), rep( 1,3), 2))

if(p == 1) {
  
  num.hessian <- hessian(fn.o, opt.o$par[-10])
  rownames(num.hessian) <- colnames(num.hessian) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3")
 
} else {
  
  num.hessian <- hessian(fn.o, opt.o$par)
  rownames(num.hessian) <- colnames(num.hessian) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b")

  }


var.ml <- solve(num.hessian)
var.matrix <-  var.ml[c(1,2,4,5,7), c(1,2,4,5,7)]

sauc.ci <- SAUC.ci(opt.o$par[c(1,2,4,5,7)], var.matrix = var.matrix, sauc.type = sauc.type)

if(p != 1) {
  
b   <- opt.o$par[10]
t33 <- (opt.o$par[6])^2
u3  <- opt.o$par[3]

f.b <- function(a){
  
  sq <- suppressWarnings(sqrt(1 + b^2 * (1 + t33 / v3)))
  
  pnorm( (a + b * u3/sqrt(v3)) / sq )
  
}


a.p <- function(a) {mean(1/f.b(a), na.rm = TRUE) - 1/p}

a.opt.try <- suppressWarnings(try(uniroot(a.p, a.interval, extendInt="yes"), silent = TRUE)) 

if (!inherits(a.opt.try, "try-error")) a.opt <- a.opt.try$root else a.opt.try <- NA

p.hat <- mean(pnorm( a.opt + b * (y3/sqrt(v3)) ), na.rm = TRUE)
p.hat2 <- 1/mean(1/f.b(a.opt), na.rm = TRUE)

}  

par <- c(opt.o$par,  a.opt, p.hat, p.hat2, sauc.ci, plogis(opt.o$par[1:2]), opt.o$convergence)
names(par) <- c( "u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3",
                 "b", "a", "p.hat","p.hat2",
                 "sauc", "sauc.lb", "sauc.ub",
                 "sen", "spe", "conv")
return(par)

}