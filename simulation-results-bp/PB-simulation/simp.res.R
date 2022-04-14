


## CR
sum(is.na(DATA[13,names(DATA[13,]) %in% "bnm.par.p"]))
sum(is.na(DATA[13,names(DATA[13,]) %in% "bnm.par.o"]))
sum(is.na(DATA[13,names(DATA[13,]) %in% "tnm.par.sa1"]))
sum(is.na(DATA[13,names(DATA[13,]) %in% "tnm.par.sa2"]))

## SAUC
summary(DATA[11,names(DATA[11,]) %in% "bnm.par.p"])
summary(DATA[11,names(DATA[11,]) %in% "bnm.par.o"])
summary(DATA[11,names(DATA[11,]) %in% "tnm.par.sa1"])
summary(DATA[11,names(DATA[11,]) %in% "tnm.par.sa2"])


summary(DATA[10,names(DATA[10,]) %in% "tnm.par.sa1"])
summary(DATA[10,names(DATA[10,]) %in% "tnm.par.sa2"])


set.seed(123)                 
aaa <- matrix(rnorm(20), 5, 4)
bbb <- cbind(aaa, aaa[, 4])   
solve(bbb)  
