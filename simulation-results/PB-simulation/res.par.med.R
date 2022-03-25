## RES.PAR


{
	conv <- matrix(DATA[13,], ncol = 8, byrow = TRUE)
	CR <- sapply(1:8, function(i){ sum(conv[,i] %in% 0)/1000 })
	dim(DATA) <- c(13,8,1000)
	
	x <- cbind(
  	apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))
	# (apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))[,1])
	)

	x[13, ] <- CR
	x <- round(x[, c(1,4,7,8)],3)
}

rownames(x) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b","sauc", "S/N/p", "CR")
colnames(x) <- c("BNM.P", "BNM.O", "SA1", "SA2")


# 
# dim(DATA) <- c(17,5,1000)
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,1]
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,2]
# 
# dim(DATA) <- c(17,5000)
