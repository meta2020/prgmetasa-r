plot(density(rnorm(1000, 0.5, 0.2)))
n.i <- 1000

x1.u <- 0.7; x1.s <- 0.1
x2.u <- 0.3; x2.s <- 0.3 
v.sd <- 0.1

X1 <- x1.u + x1.s*rlogis(n.i)
X2 <- x2.u + x2.s*rlogis(n.i)

# X <- ifelse(T.fail <= 2, X1, X2)

v.i <- rnorm(n.i, mean = 0.5, sd = v.sd)

plot(density(X1), xlim = c(-2,3))
lines(density(X2))
lines(density(v.i))


x1.u <- 0.6; x1.s <- 0.2
x2.u <- 0.4; x2.s <- 0.4
v.sd <- 0.1

X1 <- x1.u + x1.s*rlogis(n.i)
X2 <- x2.u + x2.s*rlogis(n.i)

v.i <- rnorm(n.i, mean = 0.5, sd = v.sd)

plot(density(X1), xlim = c(-2,3))
lines(density(X2))
lines(density(v.i))
