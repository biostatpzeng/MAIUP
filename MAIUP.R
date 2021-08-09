
## MAIUP depends on the HDMT function adjuated from the R package HDMT
## p1: P value for the first trait
## p2: P value for the second trait
## decorrelation=TRUE: use the decorrelation method to generate independent P vlaues for further pleiotropy analysis
## adjust=="FWER": conduct the FWER adjustment for the P_max, the corresponding output is a zero-one index vector indicating if P_max is less than the given cutoff
## adjust=="FDR": conduct the FDR adjustment for the P_max, the corresponding output is a FDR vector

##Note:
## 1. MAIUP depends on the HDMT function adjuated from the R package HDMT
## 2. MAIUP will implement quality control on the two sets of P values to remove duplicated ones; thus, to be correct, it is suggested to perform quality control before the MAIUP analysis
## 3. if any P value is zero, it will be given a value between 1E-10 and 1E-12
## 4. if any estimated proportion parameter is zero, it will be given a value of 1e-5
## 5. if any estimated proportion parameter is one,  it will be given a value of 0.9999


library(fdrtool)
MAIUP <- function(p1, p2, decorrelation=TRUE,adjust=="FWER",cutoff) {
#source("HDMT.R")
ax = cbind(p1,p2)
ax = na.omit(ax)
colnames(ax)=c("P.x","P.y")
index1 = which(ax$P.x==0)
index2 = which(ax$P.y==0)
if (length(index1)>0) {ax$P.x[index1]=runif(length(index1),1E-10,1E-12)}
if (length(index2)>0) {ax$P.x[index2]=runif(length(index2),1E-10,1E-12)}
ax = ax[!duplicated(ax$P.x),]
ax = ax[!duplicated(ax$P.y),]
px = cbind(ax$P.x,ax$P.y)

if (decorrelation==TRUE)
{
zd = cbind(qnorm(1-px[,1]),qnorm(1-px[,2]))
R <- cor_pearson(zd, px, p.threshold=1e-4)
"%^%" <- function(x, pow) {with(eigen(x), vectors %*% (values^pow * t(vectors)))}
zd <- zd %*% (R %^% (-0.5))
px <- (1-pnorm(abs(zd)))*2
}

IUT = apply(px,1,max)
nullprop <- null_estimation(input_pvalues,lambda=0.5)
if (nullprop$alpha10==0) {nullprop$alpha10=1e-5}
if (nullprop$alpha01==0) {nullprop$alpha01=1e-5}
if (nullprop$alpha00==0) {nullprop$alpha00=1e-5}
if (nullprop$alpha1==1)  {nullprop$alpha1=0.9999}
if (nullprop$alpha2==1)  {nullprop$alpha2=0.9999}

if (adjust=="FWER")
{
pnull <- adjust_quantile(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=1)
index_pleiotropy < which(IUT < quantile(pnull,probs=cutoff)) * 1 
}

if (adjust=="FDR")
{
index_pleiotropy <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,nullprop$alpha2,px,exact=1)
}

return(cbind(px,P_max=IUT,p_maiup=index_pleiotropy))
}



################### code from HDMT
adjust_quantile <-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues,exact=0){
    ## This function computes the expected quantiles of the mixture null distribution 
    ## input_pvalues is the 2-column matrix storing the two sets of p-values
    ## alpha00, alpha01, alpha10 are the estimated proportions of three nulls, 
    if (is.null(ncol(input_pvalues)))
      stop("input_pvalues should be a matrix or data frame")
    if (ncol(input_pvalues) !=2)
      stop("inpute_pvalues should have 2 column")
    input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
    if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
      warning("input_pvalues contains NAs to be removed from analysis")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
    if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
      stop("input_pvalues doesn't have valid p-values")
    #library(fdrtool)
    nmed <- nrow(input_pvalues)  
    ## compute the quantiles using the approximation method
    if (exact==0) { 
      c <- (-(1:nmed)/nmed)
      b <- alpha10+alpha01
      a <- alpha00
      pnull <- (-b+sqrt(b^2-4*a*c))/(2*a)
       
    }
    ##  compute the quantiles using the exact method
    if (exact==1) {
      cdf12 <- input_pvalues
      orderp1 <- input_pvalues[order(input_pvalues[,1]),1]
      orderp2 <- input_pvalues[order(input_pvalues[,2]),2]
      
      xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
      yy1 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit1<- gcmlcm(xx1,yy1,type="lcm")
      xknots1 <- gfit1$x.knots[-1]
      Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
      xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
      yy2 <- c(0,seq(1,nmed,by=1)/nmed)
      gfit2<- gcmlcm(xx2,yy2,type="lcm")
      xknots2 <- gfit2$x.knots[-1]
      Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
      if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
      if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))
      gcdf1 <- orderp1
      gcdf2 <- orderp2
      orderq1 <- orderp1
      orderq2 <- orderp2
      difff <- 1
      ite <- 1
      while(abs(difff)>1e-6 & ite<10) {
        #cat(ite,"..")
        for (i in 1:length(xknots1)) {
          if (i==1) {
            gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
          } else {   
            if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
              temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
              gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
            }
          }
        }
        for (i in 1:length(xknots2)) {
          if (i==1) {
            gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
          } else {   
            if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
              temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
              gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
            } 
          }
        }
        cdf12[,1] <- ifelse(gcdf1>1,1,gcdf1)
        cdf12[,2] <- ifelse(gcdf2>1,1,gcdf2)
        c <- (-(1:nmed)/nmed)
        b <- alpha10*cdf12[,1]+alpha01*cdf12[,2]
        a <- alpha00
        pnull <- (-b+sqrt(b^2-4*a*c))/(2*a)
        difff <- max(max(orderq1-pnull),max(orderq2-pnull))
        orderq1 <- pnull
        orderq2 <- pnull
        ite <- ite+1
      }
    }
    return(pnull)
  }

correct_qqplot <- function(pmax,pnull,opt="all") {
    ## pmax is the a vector of max of p-values
    ## pnull is the quantiles corresponding to pmax 
    nmed <- length(pmax)
    pindex=1:nmed
    if (opt=="subset")
    {
      ## pick a subset of p-values to plot, avoid overcrowded q-q plots###
      mpoint <- ceiling(quantile(1:length(pmax),0.95)/200)*200
      pindex <- c(1,(1:ceiling(quantile(1:length(pmax),0.95)/200))*200,mpoint:length(pmax))
    }
    
    xmax <- max(c(-log(pnull[order(pnull,decreasing=T)],base=10)[pindex],-log10(1/nmed)))
    ymax <- max(c(-log(pmax[order(pmax,decreasing=T)],base=10)[pindex],-log(pmax[order(pmax,decreasing=T)],base=10)[pindex]))
    if (xmax>0.8*ymax & xmax<1.25*ymax)
    {
      xmax <- max(xmax,ymax)
      ymax <- max(xmax,ymax)
    }
    plot((-log(pnull[order(pnull,decreasing=T)],base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],xlab="log base 10 (expected null p-values)", ylab="log base 10 (observed p-values)",col=3,xlim=c(0,xmax),ylim=c(0,ymax))
    points((-log((nmed:1)/nmed,base=10))[pindex],(-log(pmax[order(pmax,decreasing=T)],base=10))[pindex],pch=2,col=2)
    legend(0.1,max(-log(pmax,base=10)),c("Uniform null","Mixture null"),pch=2:1,col=2:3,bty="n")
    
    abline(0,1)
  }

fdr_est <-function(alpha00,alpha01,alpha10,alpha1,alpha2,input_pvalues, exact=0){
  ## alpha10,alpha01,alpha00 are estimated three types of null proportions
  ## alpha1 is the marginal null proportion for first p-value
  ## alpha2 is the marginal null proportion for second p-value
  ## input pvalues are two columns of p-values
  ## alpha is the level of FWER to be control at   
  ## check input
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")
  pmax <- apply(input_pvalues,1,max)
  nmed <- length(pmax)
  efdr1 <- rep(0,nmed)
  if (exact==0) {
   for (i in 1:nmed) {
    fdr11 <-  (pmax[i]*alpha01)/mean(pmax<=pmax[i])
    fdr12 <-  (pmax[i]*alpha10)/mean(pmax<=pmax[i])          
    fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])   
    efdr1[i] <- fdr11+fdr12+fdr2
   }  
  }
  if (exact==1) {
    #library(fdrtool)
    nmed  <- nrow(input_pvalues)  
    cdf12 <- input_pvalues
    xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
    yy1 <- c(0,seq(1,nmed,by=1)/nmed)
    gfit1<- gcmlcm(xx1,yy1,type="lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
    xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
    yy2 <- c(0,seq(1,nmed,by=1)/nmed)
    gfit2<- gcmlcm(xx2,yy2,type="lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
    if (alpha1!=1) Fknots1 <- (Fknots1 - alpha1*xknots1)/(1-alpha1) else Fknots1 <- rep(0,length(xknots1))
    if (alpha2!=1) Fknots2 <- (Fknots2 - alpha2*xknots2)/(1-alpha2) else Fknots2 <- rep(0,length(xknots2))
    orderq1 <- pmax
    orderq2 <- pmax
    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i==1) {
        gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
      } else {   
        if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
          temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
          gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
        }
      }
    }
    for (i in 1:length(xknots2)) {
      if (i==1) {
        gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
      } else {   
        if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
          temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
          gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
        } 
      }
    }
    gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
    gcdf2 <- ifelse(gcdf2>1,1,gcdf2)
    cdf12[,1] <- gcdf1
    cdf12[,2] <- gcdf2
    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*cdf12[i,2]*alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*cdf12[i,1]*alpha10)/mean(pmax<=pmax[i])          
      fdr2  <-  (pmax[i]*pmax[i]*alpha00)/mean(pmax<=pmax[i])   
      efdr1[i] <- fdr11+fdr12+fdr2
    }  
  }
  efdr1.order <- efdr1[order(pmax,decreasing=T)]
  for (i in 2:nmed)  {
    efdr1.order[i] <- min(efdr1.order[i],efdr1.order[i-1])
  }
  efdr1 <- efdr1.order[rank(-pmax)]
  return(efdr=efdr1)
}

null_estimation <- function(input_pvalues,lambda=0.5) {
  ## input_pvalues is a matrix with 2 columns of p-values, the first column is p-value for exposure-mediator association, the second column is p-value for mediator-outcome association adjusted for exposure
  ## lambda is the threshold for pi_{00} estimation, default 0.5
  #check input
  if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) !=2)
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
    stop("input_pvalues doesn't have valid p-values")
  pcut <- seq(0.1,0.8,0.1) 
  frac1 <- rep(0,8)
  frac2 <- rep(0,8)
  frac12<- rep(0,8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[,1]>=pcut[i])/(1-pcut[i])
    frac2[i] <- mean(input_pvalues[,2]>=pcut[i])/(1-pcut[i]) 
    frac12[i]<- mean(input_pvalues[,2]>=pcut[i] & input_pvalues[,1]>=pcut[i])/(1-pcut[i])^2
  }  
  ## use the median estimates for pi00 ##
  alpha00 <- min(frac12[pcut==lambda],1)
  ## alpha1 is the proportion of nulls for first p-value 
  ## alpha2 is the proportion of nulls for second p-value 
  if (ks.test(input_pvalues[,1],"punif",0,1,alternative="greater")$p>0.05) alpha1 <- 1 else   alpha1 <- min(frac1[pcut==lambda],1)  
  if (ks.test(input_pvalues[,2],"punif",0,1,alternative="greater")$p>0.05) alpha2 <- 1 else   alpha2 <- min(frac2[pcut==lambda],1)
  if (alpha00==1) {
    alpha01 <- 0
    alpha10 <- 0
    alpha11 <- 0
  } else {    
    if (alpha1==1  & alpha2==1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
      alpha00 <- 1
    }  
    if (alpha1==1  & alpha2!=1) {
      alpha10 <- 0
      alpha11 <- 0
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      alpha00 <- 1-alpha01
    }  
    if (alpha1!=1  & alpha2==1) {
      alpha01 <- 0
      alpha11 <- 0
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha00 <- 1-alpha10
    }  
    if (alpha1!=1  & alpha2!=1) {
      alpha10 <- alpha2-alpha00
      alpha10 <- max(0,alpha10)
      alpha01 <- alpha1-alpha00
      alpha01 <- max(0,alpha01)
      if ((1-alpha00-alpha01-alpha10)<0) {
        alpha11 <- 0
        alpha10 <- 1- alpha1
        alpha01 <- 1- alpha2
        alpha00 <- 1- alpha10 - alpha01
      }  else {
        alpha11 <-  1-alpha00-alpha01-alpha10
      }  
    }  
  }
  alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
  return(alpha.null)
}

fwer_est=function (alpha10, alpha01, alpha00, alpha1, alpha2, input_pvalues, 
    alpha = 0.05, exact = 0) 
{
    if (is.null(ncol(input_pvalues))) 
        stop("input_pvalues should be a matrix or data frame")
    if (ncol(input_pvalues) != 2) 
        stop("inpute_pvalues should have 2 column")
    input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
    if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues)) 
        warning("input_pvalues contains NAs to be removed from analysis")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues), 
        ]
    if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) < 
        1) 
        stop("input_pvalues doesn't have valid p-values")
    pmax <- apply(input_pvalues, 1, max)
    nmed <- length(pmax)
    c <- (-alpha)/nmed
    b <- alpha01 + alpha10
    a <- alpha00
    if (exact == 0) {
        fwer_alpha <- (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
    }
    if (exact == 1) {
        xx1 <- c(0, input_pvalues[order(input_pvalues[, 1]), 
            1])
        yy1 <- c(0, seq(1, nmed, by = 1)/nmed)
        gfit1 <- gcmlcm(xx1, yy1, type = "lcm")
        xknots1 <- gfit1$x.knots[-1]
        Fknots1 <- cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
        xx2 <- c(0, input_pvalues[order(input_pvalues[, 2]), 
            2])
        yy2 <- c(0, seq(1, nmed, by = 1)/nmed)
        gfit2 <- gcmlcm(xx2, yy2, type = "lcm")
        xknots2 <- gfit2$x.knots[-1]
        Fknots2 <- cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
        if (alpha1 != 1) 
            Fknots1 <- (Fknots1 - alpha1 * xknots1)/(1 - alpha1)
        else Fknots1 <- rep(0, length(xknots1))
        if (alpha2 != 1) 
            Fknots2 <- (Fknots2 - alpha2 * xknots2)/(1 - alpha2)
        else Fknots2 <- rep(0, length(xknots2))
        fwer_alpha <- (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
        qfwer <- fwer_alpha
        ite <- 1
        difff <- 1
        while (abs(difff) > 1e-06 & ite < 10) {
            cat(ite, "..")
            if (sum(input_pvalues[, 1] < fwer_alpha) < 60 & sum(input_pvalues[, 
                1] < fwer_alpha) > 15) {
                if (alpha1 == 1) 
                  cdf1 <- 0
                else {
                  cdf1 <- max(0, (mean(input_pvalues[, 1] < fwer_alpha) - 
                    alpha1 * fwer_alpha)/(1 - alpha1))
                  cdf1 <- min(cdf1, 1)
                }
            }
            else {
                if (sum(input_pvalues[, 1] < fwer_alpha) <= 15) 
                  cdf1 <- 1
                if (sum(input_pvalues[, 1] < fwer_alpha) >= 60) {
                  if (fwer_alpha <= xknots1[1]) 
                    cdf1 <- Fknots1[1]
                  else {
                    for (i in 2:length(xknots1)) {
                      if (sum(fwer_alpha > xknots1[i - 1] & fwer_alpha <= 
                        xknots1[i]) > 0) {
                        cdf1 <- Fknots1[i - 1] + (Fknots1[i] - 
                          Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                          1]) * (fwer_alpha - xknots1[i - 1])
                      }
                    }
                    if (fwer_alpha > xknots1[length(xknots1)]) 
                      cdf1 <- 1
                  }
                }
            }
            if (sum(input_pvalues[, 2] < fwer_alpha) < 60 & sum(input_pvalues[, 
                2] < fwer_alpha) > 15) {
                if (alpha2 == 1) 
                  cdf2 <- 0
                else {
                  cdf2 <- max(0, (mean(input_pvalues[, 2] < fwer_alpha) - 
                    alpha2 * fwer_alpha)/(1 - alpha2))
                  cdf2 <- min(cdf2, 1)
                }
            }
            else {
                if (sum(input_pvalues[, 2] < fwer_alpha) <= 15) 
                  cdf2 <- 1
                if (sum(input_pvalues[, 2] < fwer_alpha) >= 60) {
                  if (fwer_alpha <= xknots2[1]) 
                    cdf2 <- Fknots2[1]
                  else {
                    for (i in 2:length(xknots2)) {
                      if (sum(fwer_alpha > xknots2[i - 1] & fwer_alpha <= 
                        xknots2[i]) > 0) {
                        cdf2 <- Fknots2[i - 1] + (Fknots2[i] - 
                          Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                          1]) * (fwer_alpha - xknots2[i - 1])
                      }
                    }
                    if (fwer_alpha > xknots2[length(xknots2)]) 
                      cdf2 <- 1
                  }
                }
            }
            c <- (-alpha/nmed)
            if (cdf1 > 1) 
                cdf1 <- 1
            if (cdf2 > 1) 
                cdf2 <- 1
            b <- alpha10 * cdf1 + alpha01 * cdf2
            a <- alpha00
            fwer_alpha <- (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
            difff <- max(qfwer - fwer_alpha)
            qfwer <- fwer_alpha
            ite <- ite + 1
        }
    }
    return(fwer_alpha)
}
