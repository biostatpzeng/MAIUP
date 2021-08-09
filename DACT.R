#' DACT is a function used for testing for mediation effect in genome-wide epigenetic studies.
#' @param p_a, the p-value vector for the exposure-mediator associations
#' @param p_b, the p-value vector for the mediator-outcome associations
#' @param correction, the name of method for correction, either "Efron" or "JC". By default, it is set to NULL.
#' @return the p-value vector of mediation effect testing
#' @references Liu, Z. , Shen, J., Barfield, R., Schwartz, J., Baccarelli, A., Lin, X., 2020. Large-scale hypothesis testing for causal mediation effects in Genome-wide epigenetic studies. (under review)
#' @author Zhonghua Liu < zhhliu@hku.hk >
#' @usage DACT(p_a,p_b,correction)
#' @export
DACT = function(p_a,p_b,correction="NO"){
  Z_a = stats::qnorm(p_a,lower.tail = F)
  Z_b = stats::qnorm(p_b,lower.tail = F)
  pi0a = nonnullPropEst(Z_a,0,1)
  pi0b = nonnullPropEst(Z_b,0,1)
  #pi0a = locfdr::locfdr(Z_a,nulltype = 0)$fp0[5,3]
  #pi0b = locfdr::locfdr(Z_b,nulltype = 0)$fp0[5,3]
  if(pi0a > 1){
    pi0a = 1
  }
  if(pi0b >1){
    pi0b = 1
  }
  p.mat = cbind(p_a,p_b)
  p3 = (apply(p.mat,1,max))^2
  wg1 = pi0a*(1-pi0b)
  wg2 = (1-pi0a)*pi0b
  wg3 = pi0a*pi0b
  wg.sum = wg1 + wg2 + wg3
  wg.std = c(wg1,wg2,wg3)/wg.sum
  p_dact = wg.std[1]*p_a + wg.std[2]*p_b + wg.std[3]*p3
  if (correction == "NO")  {p_dact = p_dact}
  if (correction == "Efron") {p_dact = EfronCorrect(p_dact)}
  if (correction == "JC") {p_dact = JCCorrect(p_dact)}
  return(p_dact)
}

EfronCorrect = function(pval){
  z = stats::qnorm(1-pval)
  res <- locfdr(z,nulltype = 1,plot=0)
  mean.emp = res$fp0["mlest","delta"]
  sd.emp = res$fp0["mlest","sigma"]
  pval.emp = stats::pnorm(z,mean = mean.emp,sd = sd.emp,lower.tail = F)
  return(pval.emp)
}

JCCorrect = function(pval){
  z = stats::qnorm(pval,lower.tail = F)
  res= nullParaEst(z)
  pval.JC = stats::pnorm(z,mean = res$mu,sd = res$s,lower.tail = F)
  return(pval.JC)
}

nonnullPropEst <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  epsest=NULL
  for (j in 1:length(tt)) {
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    }
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

nullParaEst<-function (x,gamma=0.1)
{
 # x is a vector of z-values
 # gamma is a parameter, default is 0.1
 # output the estimated mean and standard deviation
 n = length(x)
 t = c(1:1000)/200
 gan    = n^(-gamma)
 that   = 0
 shat   = 0
 uhat   = 0
 epshat = 0
 phiplus   = rep(1,1000)
 phiminus  = rep(1,1000)
 dphiplus  = rep(1,1000)
 dphiminus = rep(1,1000)
 phi       = rep(1,1000)
 dphi      = rep(1,1000)
 for (i in 1:1000) {
    s = t[i]
    phiplus[i]   = mean(cos(s*x))
    phiminus[i]  = mean(sin(s*x))
    dphiplus[i]  = -mean(x*sin(s*x))
    dphiminus[i] = mean(x*cos(s*x))
    phi[i]       = sqrt(phiplus[i]^2 + phiminus[i]^2)
 }
 ind = min(c(1:1000)[(phi - gan) <= 0])
 tt = t[ind]
 a  = phiplus[ind]
 b  = phiminus[ind]
 da = dphiplus[ind]
 db = dphiminus[ind]
 c  = phi[ind]
 that   = tt
 shat   = -(a*da + b*db)/(tt*c*c)
 shat   = sqrt(shat)
 uhat   = -(da*b - db*a)/(c*c)
 epshat = 1 - c*exp((tt*shat)^2/2)
 return(musigma=list(mu=uhat,s=shat))
}
