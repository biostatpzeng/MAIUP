

library(MASS)
library(fdrtool) ##gcmlcm --> HDMT
library(locfdr) ##efron --> DACT
library(data.table)
source("DACT.R")
source("PLACO.R")


N = 1000
FWER = 0.01
pai = rbind(
c(0.60,0.20,0.20,0.00),
c(0.90,0.05,0.05,0.00),
c(1.00,0.00,0.00,0.00),
c(0.40,0.20,0.20,0.20),
c(0.80,0.05,0.05,0.10))
pai=pai[pai_x,]
taux = c(2,3,4)[tau_x]
S = c(10000,15000,20000)[SSS_x]
rou = 0
Sigma=as.matrix(cbind(c(1,rou),c(rou,1)))
res1 = matrix(NA, N, 6)

mu1=c(0,0)
mu2=c(taux,0)
mu3=c(0,taux)
mu4=c(taux,taux)

for (j in 1:N)
{
try({
set.seed(j)
if (pai_x<=2) {
zx=rbind(mvrnorm(S*pai[1],mu1,Sigma,empirical=TRUE),mvrnorm(S*pai[2],mu2,Sigma,empirical=TRUE),mvrnorm(S*pai[3],mu3,Sigma,empirical=TRUE))
}
if (pai_x==3) {zx=mvrnorm(S*pai[1],c(0,0),Sigma,empirical=TRUE)}
if (pai_x>3) {
zx=rbind(mvrnorm(S*pai[1],mu1,Sigma,empirical=TRUE),mvrnorm(S*pai[2],mu2,Sigma,empirical=TRUE),mvrnorm(S*pai[3],mu3,Sigma,empirical=TRUE),mvrnorm(S*pai[4],mu4,Sigma,empirical=TRUE))
}

p1 = pnorm(-abs(zx[,1]))*2
p2 = pnorm(-abs(zx[,2]))*2
input_pvalues <- cbind(p1,p2)

px=input_pvalues
IUT = apply(input_pvalues,1,max)  ###### IUT

zd=zx
index=c(which(is.infinite(zd[,1])),which(is.infinite(zd[,2])),which(is.nan(zd[,1])),which(is.nan(zd[,2])),which(px[,1]==0),which(px[,2]==0))
index=unique(index)
if (length(index)>0) {px=px[-index,];zd=zd[-index,]}

pd1 = DACT(p_a=px[,1],p_b=px[,2],correction="NO") # DACT_NO
pd2 = DACT(p_a=px[,1],p_b=px[,2],correction="JC") # DACT_JC
pd3 = DACT(p_a=px[,1],p_b=px[,2],correction="Efron") # DACT_Efron

###### placo
varz=var_placo(zd, px, p.threshold=-1)
out <- t(sapply(1:dim(zd)[1], function(i) placo(Z=zd[i,], VarZ=varz)))
p_placo = out[,2]
###### placo

p_MAIUP <- (px[,1], px[,2], decorrelation=FALSE,adjust=="FWER",FWER)###### MAIUP
res1[j,]=c(apply(cbind(IUT,p_placo,pd1,pd2,pd3)<FWER,2,mean),mean(p_MAIUP))
})
}
