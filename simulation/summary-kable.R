##
## SUMMARY
##



##
## SUMMARY OF SAUC ----
##
##
## TABLE 1 ----
##
.table.sauc <- function(Sn=1,npt.n=1,p=7){
  
  res = NULL
  for(mark.n in 1:5){
    
  load(paste0(p,"npt/simipd_npt",npt.n,"_marker",mark.n ,"_S",Sn,".RData"))
  DATA[c(1:10, 12), DATA[11,]!=0] <- NA
  dim1 <- nrow(DATA)
  dim2 <- 3
  
  name.c <- colnames(DATA)[1:dim2]
  
  dim(DATA) <- c(dim1, dim2, 1000)
  
  med <- apply(DATA, 1:2, function(x) median(x, na.rm = TRUE))[12,]
  q1  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.25, na.rm = TRUE))[12,]
  q3  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.75, na.rm = TRUE))[12,]
  conv<- 1-apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[11,3]
  sn = apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[13,]
  
  s=c(sprintf("%1.f (%1.f)",sn[1], sn[2]))
  sauc = c(sprintf("%.2f (%.2f, %.2f)",med, q1, q3))
  res = rbind(res, c(s, mark.n, sauc, round(conv*100,2), sprintf("%.2f", med[2]-med[1]), sprintf("%.2f", med[3]-med[1])))
  }
  
  
  colnames(res) <- c("$S (N)$","B","BNM$_P$", "BNM$_O$", "Proposal", "CR", "PB", "Bias")
  
  # gsub("NaN \\(NA\\)", "", res)
  
  return(res)
}

tab.npt1 <- rbind(
  .table.sauc(1,1,7),
  .table.sauc(2,1,7),
  .table.sauc(3,1,7))

tab.npt2 <- rbind(
  .table.sauc(1,2,7),
  .table.sauc(2,2,7),
  .table.sauc(3,2,7))


tab.npt1 <- rbind(
  .table.sauc(4,1,5),
  .table.sauc(5,1,5),
  .table.sauc(6,1,5))

tab.npt2 <- rbind(
  .table.sauc(4,2,5),
  .table.sauc(5,2,5),
  .table.sauc(6,2,5))


npt = c("50-150", rep("", 14), "40-300", rep("", 14))
sauc.all = rbind(tab.npt1, tab.npt2)
sauc.all[-seq(1,30,5),1] <- ""
sumtab = cbind(npt,sauc.all)
library(kableExtra)
sumtab[,1:7] %>%
  kbl(
    caption = "Comparison of medians with the first and third quantiles of estimates of SAUC(2) by the HZ model and the proposed method.",
    # format = ifelse(save.tex, "latex", "html"),
    format = "latex",
    longtable = F, 
    booktabs = T,
    align = "r",
    position = "!htb",
    linesep = c("", "", "", "", "\\addlinespace"),
    escape = FALSE,
    label = "tab:sum1")%>%
  footnote(general = "
			B denotes the scenarios of biomarker coresponding to Table 2;
			CR shows convergence rate of the proposed method;
			estimates are summarized by median (first quantile, third quantiles).", 
           escape = FALSE, 
           threeparttable = TRUE,  
           general_title = "") 

##
## TABLE 2 ----
##
.table.sauc.mis <- function(Sn,npt.n, censd){
  res = NULL
  for(mark.n in 1:5){
    
    load(paste0("nptmis-", censd, "/simipd_npt",npt.n,"_marker",mark.n ,"_S",Sn,".RData"))
    DATA[c(1:10, 12), DATA[11,]!=0] <- NA
    dim1 <- nrow(DATA)
    dim2 <- 3
    
    name.c <- colnames(DATA)[1:dim2]
    
    dim(DATA) <- c(dim1, dim2, 1000)
    
    med <- apply(DATA, 1:2, function(x) median(x, na.rm = TRUE))[12,]
    q1  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.25, na.rm = TRUE))[12,]
    q3  <- apply(DATA, 1:2, function(x) quantile(x, probs = 0.75, na.rm = TRUE))[12,]
    conv<- 1-apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[11,3]
    sn = apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))[13,]
    
    s=c(sprintf("%1.f (%1.f)",sn[1], sn[2]))
    sauc = c(sprintf("%.2f (%.2f, %.2f)",med, q1, q3))
    res = rbind(res, c(s, mark.n, sauc, round(conv*100,2), sprintf("%.2f", med[2]-med[1]), sprintf("%.2f", med[3]-med[1])))
  }
  
  
  colnames(res) <- c("S (N)","Biomarker","BNM_P", "BNM_O", "Proposal", "CR", "PB", "Bias")
  
  # gsub("NaN \\(NA\\)", "", res)
  
  return(res)
}

tab.npt1 <- rbind(
  .table.sauc.mis(1,1,"LN"),
  .table.sauc.mis(2,1,"LN"),
  .table.sauc.mis(3,1,"LN"))
tab.npt2 <- rbind(
  .table.sauc.mis(1,2,"LN"),
  .table.sauc.mis(2,2,"LN"),
  .table.sauc.mis(3,2,"LN"))

tab.npt1 <- rbind(
  .table.sauc.mis(1,1,"UNIF"),
  .table.sauc.mis(2,1,"UNIF"),
  .table.sauc.mis(3,1,"UNIF"))
tab.npt2 <- rbind(
  .table.sauc.mis(1,2,"UNIF"),
  .table.sauc.mis(2,2,"UNIF"),
  .table.sauc.mis(3,2,"UNIF"))

npt = c("50-150", rep("", 14), "40-300", rep("", 14))
sauc.all = rbind(tab.npt1, tab.npt2)
sauc.all[-seq(1,30,5),1] <- ""
sumtab = cbind(npt,sauc.all)
library(kableExtra)
sumtab[,1:7] %>%
  kbl(
    caption = "Comparison of medians with the first and third quantiles of estimates of SAUC(2) by the HZ model and the proposed method.",
    # format = ifelse(save.tex, "latex", "html"),
    format = "latex",
    longtable = F, 
    booktabs = T,
    align = "r",
    position = "!htb",
    linesep = "",
    escape = FALSE,
    label = "tab:sum1")%>%
  footnote(general = "
			B denotes the scenarios of biomarker coresponding to Table 2;
			CR shows convergence rate of the proposed method;
			estimates are summarized by median (first quantile, third quantiles).", 
           escape = FALSE, 
           threeparttable = TRUE,  
           general_title = "") 


.table.pars <- function(Sn, mark.n){
  
  load(paste0("npt1/simipd_npt1_marker",mark.n ,"_S",Sn,".RData"))
  # DATA[c(1:10, 12), DATA[11,]!=0] <- NA
  dim1 <- nrow(DATA)
  dim2 <- ncol(DATA)/1000
  
  name.c <- colnames(DATA)[1:dim2]
  
  dim(DATA) <- c(dim1, dim2, 1000)
  
  med <- apply(DATA, 1:2, function(x) mean(x, na.rm = TRUE))
  sd  <- apply(DATA, 1:2, function(x) sd(x, na.rm = TRUE))
  rownames(med) <- c("mu1","mu2","mu3", "tau1","tau2","tau3", "rho1","rho2","rho3", "beta","conv", "sauc","SNp")
  
  ## conv
  med[11,] = 1-med[11,]
  s=c(sprintf("%1.f (%1.f)",med[13,1], med[13,2]), rep("",nrow(med)-1))
  res = cbind(s, rownames(med), matrix(sprintf("%.3f (%.3f)",med, sd), ncol=3))
  
  colnames(res) <- c("S (N)","Par","BNM_P", "BNM_O", "Proposal")
  
  gsub("NaN \\(NA\\)", "", res)
  
  return(res[1:10,])
}


rbind(
  .table.pars(1,1),
  .table.pars(2,1),
  .table.pars(3,1))


##
## ALPHA
##

load("5npt/alpha-npt1.RData")
alpha1 = alpha
load("5npt/alpha-npt2.RData")
alpha2 = alpha
load("7npt/alpha-npt1.RData")
alpha11 = alpha
load("7npt/alpha-npt2.RData")
alpha21 = alpha
load("nptmis-UNIF/alpha-npt1.RData")
alpha1mis1 = alpha
load("nptmis-UNIF/alpha-npt2.RData")
alpha2mis1 = alpha
load("nptmis-LN/alpha-npt1.RData")
alpha1mis2 = alpha
load("nptmis-LN/alpha-npt2.RData")
alpha2mis2 = alpha

atab = round(rbind.data.frame(alpha1, alpha2, alpha11, alpha21, alpha1mis1, alpha2mis1, alpha1mis2, alpha2mis2),2)
colnames(atab) = paste0("Biomarker",1:5)
atab$Patients = rep(c("50-150", "40-300"),4)
atab$Consor = c("Exp(0.2)","","","","U[1,4]","","LN","")
atab$p = c("0.7","","0.5","","0.7","","0.7","")
atab[,c(7,8,6,1:5)]%>%
  kbl(
    caption = "The values of $\\alpha$ in different simulation scenarios.",
    # format = ifelse(save.tex, "latex", "html"),
    format = "latex",
    longtable = F, 
    booktabs = T,
    align = "r",
    position = "!htb",
    linesep = "",
    escape = FALSE,
    label = "tab:sum1")%>%
  footnote(general = "
			Censor indicates the true censoring distribution.", 
           escape = FALSE, 
           threeparttable = TRUE,  
           general_title = "") 
