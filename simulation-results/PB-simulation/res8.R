## RES


{
	
	conv <- matrix(DATA[13,], ncol = 8, byrow = TRUE)
	CR <- sapply(1:8, function(i){ sum(conv[,i] %in% 0)/1000 })
	
	
	p <- matrix(DATA[12,], ncol = 8, byrow = TRUE)
	n.p <- colMeans(p, na.rm = TRUE)
	
	
	sauc <- matrix(DATA[11,], ncol = 8, byrow = TRUE)
	boxplot(sauc, outline =TRUE, col=c("lightgray", "lightgray", "white", "white", "lightgray", "lightgray"), ylim = c(0.5,1),
									names = c("BNM.P", "TNM.P", "TNM.P2", "BNM.O", "TNM.O", "TNM.O2", "SA1", "SA2"), boxwex = 0.5)
	
	
	med  <- sapply(1:8, function(i) median(sauc[,i], na.rm = TRUE))
	abline(h = med[1], col = 2, lwd=2, lty=2)
	mean <- colMeans(sauc, na.rm = TRUE)
	points(mean, type = "o", col=2, pch=19)
	# sd   <- sapply(1:4, function(i) sd(sauc[,i], na.rm = TRUE))
	
	res <- med 
}

x <- round(rbind(res, res - res[c(1,2,3,1,2,3,1,1)], CR, n.p),3)
# x[3:4, c(1,4)] <- NA; x[5:6, c(1,2,3)] <- NA
rownames(x) <- c("sauc.Med", "Bias.Med", "CR", "S/N/p")
colnames(x) <- c("BNM.P", "TNM.P", "TNM.P2", "BNM.O", "TNM.O", "TNM.O2", "SA1", "SA2")

# 
# dim(DATA) <- c(17,5,1000)
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,1]
# apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))- apply(DATA, 1:2,function(x) mean(x, na.rm=TRUE))[,2]
# 
# dim(DATA) <- c(17,5000)
