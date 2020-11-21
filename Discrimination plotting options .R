
library(survival)
library(survC1)
library(survIDINRI)
data(pbc)

pbc <- within(pbc, {
  ## transplant (1) and death (2) are considered events, and marked 1
  event <- as.numeric(status %in% c(1,2))

  ## Create a survival vector
  Surv <- Surv(time, event)
})

## Create numeric variable for sex
pbc$female <- as.numeric(pbc$sex == "f")

## C-statistic for age sex model
unoC.age.female <- Est.Cval(mydata = pbc[,c("time","event","age","female")],tau = 10 * 365.25)
unoC.age.female$Dhat
unoC.age.female.albumin <- Est.Cval(mydata = pbc[,c("time","event","age","female","albumin")], tau = 10 * 365.25)
unoC.age.female.albumin$Dhat
uno.C.delta <- Inf.Cval.Delta(mydata = pbc[,c("time","event")],
                              covs0  = pbc[,c("age","female")],           # age sex model
                              covs1  = pbc[,c("age","female","albumin")], # age sex albumin model
                              tau    = 10 * 365.25,                       # Trucation time (max time to consider)
                              itr = 1000)                                 # Iteration of perturbation-resampling
uno.C.delta


res.IDI.INF <- IDI.INF(indata = pbc[,c("time","event")],
                       covs0 = pbc[,c("age","female")],           # age sex model
                       covs1 = pbc[,c("age","female","albumin")], # age sex albumin model
                       t0 = 10 * 365.25,
                       npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)

## M1 IDI; M2 continuous NRI; M3 median improvement
IDI.INF.OUT(res.IDI.INF)

## M1 red area; M2 distance between black points; M3 distance between gray points
IDI.INF.GRAPH(res.IDI.INF)

# package survivalROC

install.packages("survivalROC")
library(survivalROC)
data(mayo)
nobs <- NROW(mayo)
cutoff <- 365
Staltscore4 <- NULL
Mayo.fit4 <- survivalROC.C( Stime = mayo$time,
                            status = mayo$censor, marker = mayo$mayoscore4, predict.time = cutoff, span = 0.25*nobs^(-0.20))
Staltscore4 <- Mayo.fit4$Survival
plot(Mayo.fit4$FP, Mayo.fit4$TP, type = "l",
     xlim = c(0,1), ylim = c(0,1),
     xlab = paste( "FP \n AUC =",round(Mayo.fit4$AUC,3)), ylab = "TP",main = "Year = 1" )
abline(0,1)
