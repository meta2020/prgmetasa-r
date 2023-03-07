##
## 2. SELECT INDICATOR (TRUE OR FALSE)
##

select.id <- function(p_dataSt.df, b, a) {
	
	select.p  <- pnorm(b * p_dataSt.df$t_lnHR + a)
	select.id <- sapply(1:nrow(p_dataSt.df), function(i){rbinom(1, 1, select.p[i]) })
	
	# s.id <- (select.id==1)
	
	select.id
}
