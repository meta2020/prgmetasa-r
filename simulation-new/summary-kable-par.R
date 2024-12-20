##
## SUMMARY
##

library(kableExtra)

##
## SUMMARY OF SAUC ----
##
##
## TABLE 1 ----
##
.table.par <- function(Sn=1,npt.n=1,p=7){  ## res.row = 1(logit-sensitivity);res.row = 2(logit-specificity)
  
  res = NULL
  for(mark.n in 1:6){
    
    load(paste0(p,"npt/simipd_npt",npt.n,"_marker",mark.n ,"_S",Sn,".RData"))
    DATA[c(1:10, 1:2), DATA[11,]!=0] <- NA
    dim1 <- nrow(DATA)
    dim2 <- 5
    
    name.c <- colnames(DATA)[1:dim2]
    
    dim(DATA) <- c(dim1, dim2, 1000)
    
    med <- apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[1:9,]
    q1  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.25, na.rm = TRUE))[1:9,]
    q3  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.75, na.rm = TRUE))[1:9,]
    # conv<- 1-apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[11,3]
    sn = apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[13,]
    
    s=sprintf("%1.f (%1.f)",sn[1], sn[3])
    sauc=matrix(sprintf("%.2f (%.2f, %.2f)",med, q1, q3), nrow=9, ncol=5)[1:2,]
    res = rbind(res, cbind(rbind(s,""), c(mark.n,""), c("$\\mu_{se}$","$\\mu_{sp}$"),sauc,
                           sprintf("%.2f", (med[1:2,3]-med[1:2,1])), sprintf("%.2f", (med[1:2,5]-med[1:2,1]))
    ))
    # sprintf("%.2f", 100*(med[3]-med[2])/(med[2]-med[1]))))
  }
  
  
  colnames(res) <- c("$S (N)$","B","Par","BNM$_P$", "TNM$_P$", "BNM$_O$", "TNM$_O$", "Proposal","RB","Bias")
  rownames(res) = NULL
  # gsub("NaN \\(NA\\)", "", res)
  
  return(res)
}

tab.par1 = rbind(
  .table.par(1,1,7),
  .table.par(2,1,7),
  .table.par(3,1,7))

tab.par2 = rbind(
  .table.par(1,2,7),
  .table.par(2,2,7),
  .table.par(3,2,7))


tab.par1 <- rbind(
  .table.par(4,1,5),
  .table.par(5,1,5),
  .table.par(6,1,5))

tab.par2 <- rbind(
  .table.par(4,2,5),
  .table.par(5,2,5),
  .table.par(6,2,5))


npt = c("50-150", rep("", 35), "40-300", rep("", 35))
sauc.all = rbind(tab.par1, tab.par2)
sauc.all[-seq(1,72,12),1] <- ""
sumtab = cbind(npt,sauc.all)

# sumtab[1:36,c(1:5,7,9:11)] %>%
sumtab[37:72,c(1:5,7,9:11)] %>%
  kbl(
    caption = "Comparison of medians with the first and third quantiles of estimates of SAUC(2) by the HZ model and the proposed method.",
    # format = ifelse(save.tex, "latex", "html"),
    format = "latex",
    longtable = F, 
    booktabs = T,
    align = "r",
    position = "!htb",
    linesep = c(rep("",11), "\\addlinespace"),
    escape = FALSE,
    label = "tab:sum1")%>%
  footnote(general = "
			B denotes the scenarios of biomarker coresponding to Table 2;
			CR shows convergence rate of the proposed method;
			estimates are summarized by median (first quantile, third quantiles);
			RB denotes reporting bias;
      Bias denotes bias of the proposed method.", 
           escape = FALSE, 
           threeparttable = TRUE,  
           general_title = "") 

