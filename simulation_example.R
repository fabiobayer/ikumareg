# Author: Fabio M Bayer
# Data: March 08, 2021
# e-mail: bayer@ufsm.br

source("ikr.fit.R") # read the function which fits the inflated Kumaraswamy regression model

R = 10000 # Monte Carlo replications
vn = c(1000,500,200,100,50) # sample sizes
alfa = 0.05 #  coverage rate (1-alfa)
set.seed(2021)

# parameter values
vgamma = c(-0.5,-1) 
vpi = c(1,-1)
vbeta = c(1,-2)
vsigma = c(1,1.5)

pa = c(vgamma,vpi,vbeta,vsigma)

for(n in vn)
{
  # covariates 
  z = as.matrix(runif(n))
  w = as.matrix(runif(n))
  x = as.matrix(runif(n))
  q = as.matrix(runif(n))

  # with intercepts
  Z = cbind(1,z)
  W = cbind(1,w)
  X = cbind(1,x)
  Q = cbind(1,q)
  
  # linear predictors
  eta1 = Z %*% vgamma
  eta2 = W %*% vpi
  eta3 = X %*% vbeta
  eta4 = Q %*% vsigma
  
  lambda = linkinv1(eta1)
  p = linkinv2(eta2)
  omega = linkinv3(eta3)
  phi = linkinv4(eta4)
  
  # inicializations
  TC<-c() 
  coef_results<-c()
  results<-c()
  
  nc=0 # nonconvergence
  ni=0 # problems with information matrix inversion
  
  i = 0
  while(i<R) # Monte Carlo loop
  {
    y = rikum(n,lambda = lambda, p=p, omega = omega, phi = phi) # random generation
    fit <- ikr.fit(z,w,x,q,y,diag=0) # fit the model
    
    if( (fit$conv==0)  && (fit$cov_ok==0) ) 
    {
      results<-rbind(results,fit$pvalues)
      coef_results<-rbind(coef_results,fit$coef)
      
      LI<- fit$coef - qnorm(1-alfa/2)* fit$stderror
      LS<- fit$coef + qnorm(1-alfa/2)* fit$stderror
      
      TC <- rbind(TC, ((pa<LI) + (pa>LS)) )
      
      i = i+1
    }else{ # nonconvergence failures
      if(fit$cov_ok==1) ni = ni+1
      nc<-nc+1
      print(c("Nonconvergence",i,nc))
    }
  }
  
  m<-colMeans((coef_results),na.rm=T)
  sd<-apply(coef_results,2,sd)
  bias<-m-pa
  rb<- 100*bias/pa
  tc <- 1-colMeans((TC),na.rm=T) #  coverage rate
  
  M<-rbind(pa,m,sd,bias,rb,tc)
  row.names(M)<-c("Parameters","Mean","EP","Bias","RB","CR")
  
  print(c("n",n),quote=F)
  print(c("nc",nc-ni),quote=F)
  print(c("ni",ni),quote=F)
  print(round(M,3))
}
