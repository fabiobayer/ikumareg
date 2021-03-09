# Author: Fabio M Bayer
# Data: March 08, 2021
# e-mail: bayer@ufsm.br

library(extraDistr)

# density function
dikum<-function(y,lambda=0.2,p=0.2,omega=0.5,phi=5)
{
  f0 = ifelse((y == 0), lambda*(1-p), 0)
  f1 = ifelse((y == 1), lambda*p, 0)
  
  f =  ifelse((y == 0), 0, (1-lambda)*dkumar(y,phi,log(0.5)/log(1-omega^phi)))
  f =  ifelse((y == 1), 0, f)

  (f0 + f1 + f)
}

# cumulative distribution function
pikum<-function(y,lambda=0.2,p=0.2,omega=0.5,phi=5)
{
  p0 = lambda*(1-p)
  p1 = ifelse((y == 1), lambda*(p), 0)
  p<- (1-lambda)*(1-(1-(y^phi))^( (log(0.5))/(log(1-(omega^phi))) ))
  ret = ifelse(y==1, p0+p1+p, p0+p)
  ret
}

# quantile function
qikum<-function(u,lambda=0.2, p=0.2, omega=0.5,phi=5)
{
  n = length(u)
  if(length(p)==1) p <- rep(p,n)

  if(n>1)
  {
    ret=rep(NA,n)
    for(i in 1:n)
    {
      if(u[i] < (lambda[i]*(1-p[i])) ) 
      {
        ret[i]=0
      }else{
        if(u[i] > (1-lambda[i]*p[i]) ) 
        {
          ret[i]=1
        }else{
          ret[i]= qkumar((u[i]-lambda[i]*(1-p[i]))/(1-lambda[i]), phi[i], log(0.5)/log(1-omega[i]^phi[i]))
        }
      }
    }
  }else{
    if(u < (lambda*(1-p)) ) 
    {
      ret=0
    }else{
      if(u > (1-lambda*p) ) 
      {
        ret=1
      }else{
        ret=qkumar((u-lambda*(1-p))/(1-lambda), phi, log(0.5)/log(1-omega^phi))
      }
    }
  }
  ret
}

# inversion method for randon generation
rikum<-function(n,lambda=0.2,p=0.2,omega=0.5,phi=5)
{
  u <- runif(n)
  qikum(u,lambda=lambda,p=p,omega=omega,phi=phi)
}

