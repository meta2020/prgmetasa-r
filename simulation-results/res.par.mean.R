## RES.PAR


{
	dim(DATA) <- c(13,4,1000)
	
	x <- cbind(
  	apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE)),
	(apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,1])
	)

	x <- round(x[-c(12,13), -5],3)
}

rownames(x) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b","sauc")
colnames(x) <- c("BNM.P", "BNM.O", "SA1", "SA2", "PB", "Bias1", "Bias2" )


# 
# dim(DATA) <- c(17,5,1000)
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,1]
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,2]
# 
# dim(DATA) <- c(17,5000)
