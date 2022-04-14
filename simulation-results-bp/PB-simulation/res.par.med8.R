## RES.PAR


{
	conv <- matrix(DATA[13,], ncol = 8, byrow = TRUE)
	CR <- sapply(1:8, function(i){ sum(conv[,i] %in% 0)/1000 })
	
	
	dim(DATA) <- c(13,8,1000)
	
	x <- cbind(
  	apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))
	# (apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) median(x, na.rm=TRUE))[,1])
	)

	x <- round(x,3)
	x[13, ] <- CR
}

rownames(x) <- c("u1", "u2", "u3", "t1", "t2", "t3", "r1", "r2", "r3", "b","sauc", "S/N/p", "conv")
colnames(x) <- c("BNM.P", "TNM.Cho", "TNM.exp","BNM.O", "TNM.Cho.O", "TNM.exp.O", "SA1", "SA2")


# 
# dim(DATA) <- c(17,5,1000)
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,1]
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,2]
# 
# dim(DATA) <- c(17,5000)
