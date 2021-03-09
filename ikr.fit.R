# Author: Fabio M Bayer
# Data: March 08, 2021
# e-mail: bayer@ufsm.br

source("auxiliar.R") 
source("ikum-omega-phi.R") 
library(rootSolve)

"ikr.fit" <- function(z,w,x,q,y,diag=1,theta0=NA){
  
  k <- c()
  
  k$y <- y
  k$z <- z
  k$w <- w
  k$x <- x
  k$q <- q
  
  z <- as.matrix(data.matrix(z))
  z1 <- cbind(1,z)
  w <- as.matrix(data.matrix(w))
  w1 <- cbind(1,w)
  x <- as.matrix(data.matrix(x))
  x1 <- cbind(1,x)
  q <- as.matrix(data.matrix(q))
  q1 <- cbind(1,q)
  y <- as.vector(y)

  m <- ncol(z1)
  u <- ncol(w1)
  r <- ncol(x1)
  s <- ncol(q1)
  n <- length(y)
  n0<- sum(y==0) # number of zeros
  n1<- sum(y==1) # number of ones
  
  if(n0 == 0) 
  {
    p=rep(1,n) 
    k$p <- p
    u=0
  }
    
  if(n1 == 0)
  {
    p=rep(0,n)
    k$p <- p
    u=0
  }
  
  i01  <- (y==0)+(y==1)
  
  # para teste sob parametros
  if(is.na(theta0[1]))
  {
    theta0<-rep(0,(m+u+r+s))
  }
  eta <- c()
 
  loglik <- function(theta) 
  {
    gama <- as.vector(theta[1:m])
    beta <- as.vector(theta[(m+u+1):(m+u+r)])
    sigma <- as.vector(theta[(m+u+r+1):(m+u+r+s)])

    eta1 <- as.vector(z1%*%gama)
    eta3 <- as.vector(x1%*%beta)
    eta4 <- as.vector(q1%*%sigma)
    
    lambda <- as.vector(linkinv1(eta1))
    omega <- as.vector(linkinv3(eta3))
    phi <- as.vector(linkinv4(eta4))
    
    l2 = 0
    if(u!=0)
    {
      pi <- as.vector(theta[(m+1):(m+u)])
      eta2 <- as.vector(w1%*%pi)
      p <- as.vector(linkinv2(eta2))
      
      l2a = ifelse((y == 1), log(p), 0)
      l2b = ifelse((y == 0), log(1-p), 0)
      l2 = l2a + l2b
    }

    l1 = ifelse((i01 == 1), log(lambda), log(1-lambda))
    
    l3 = ifelse((i01 == 1), 0, 
                log(phi)+
                  log(log(0.5)/log(1-omega^phi))+
                  (phi-1)*log(y)+
                  ( (log(0.5)/log(1-omega^phi))-1)*log(1-y^phi)
    )
    
    sum(l1 + l2 + l3)
  }
  
  escore <- function(theta) 
  {
    gama <- as.vector(theta[1:m])
    beta <- as.vector(theta[(m+u+1):(m+u+r)])
    sigma <- as.vector(theta[(m+u+r+1):(m+u+r+s)])
    
    eta1 <- as.vector(z1%*%gama)
    eta3 <- as.vector(x1%*%beta)
    eta4 <- as.vector(q1%*%sigma)
    
    lambda <- as.vector(linkinv1(eta1))
    omega <- as.vector(linkinv3(eta3))
    phi <- as.vector(linkinv4(eta4))
    
    if(u!=0)
    {
      pi <- as.vector(theta[(m+1):(m+u)])
      eta2 <- as.vector(w1%*%pi)
      p <- as.vector(linkinv2(eta2))
    }
    
    delta = log(0.5)/log(1-omega^phi)
    ci = ((omega^(phi-1))/((1-omega^phi)*log(1-omega^phi)))*(delta*log(1-y^phi)+1)
      
    a = ifelse( (y == 0) , (1)/(lambda), (-1)/((1-lambda)))
    a = ifelse( (y == 1) , (1)/(lambda), a)
    
    rho1 = ifelse( (y == 1), 1/p, 0 )
    rho0 = ifelse( (y == 0), -1/(1-p), 0 )
    rho = rho1 + rho0
    
    c = ifelse((y == 0), 0, phi*ci)
    c = ifelse((y == 1), 0, c)
    
    v = ifelse((y == 0), 0, 1/phi + log(y) + ci*omega*log(omega) -
                 (delta-1)*((y^phi)*log(y)/(1-y^phi)))
    v = ifelse((y == 1), 0, v)
    
    T1 = diag(as.vector( 1/diflink1(lambda) )) 
    T2 = diag(as.vector( 1/diflink2(p) )) 
    T3 = diag(as.vector( 1/diflink3(omega) )) 
    T4 = diag(as.vector( 1/diflink4(phi) )) 
    
    Ugama <- t(z1) %*% T1 %*% a 
    Upi <- t(w1) %*% T2 %*% rho 
    Ubeta <- t(x1) %*% T3 %*% c 
    Usigma <- t(q1) %*% T4 %*% v 
    
    if(u != 0) 
    {
      derivada=c(Ugama,Upi,Ubeta,Usigma)
    }else{
      derivada=c(Ugama,Ubeta,Usigma)
    }
    derivada
  }
  
  ### initial values
  ys01<- y
  if(any(y==0)) ys01<- ys01[-which(y==0)]
  if(any(y==1)) ys01<- ys01[-which(ys01==1)]
  
  y01<-y
  y01[which(y01==0)]<- min(ys01)
  y01[which(y01==1)]<- max(ys01)

  #print(x1)
  ystar = linkfun3(y01)
  beta_ini <- (lm.fit(x1, ystar))$coefficients
  #print(beta_ini)
  
  gama_ini = rep(0,m)
  sigma_ini = rep(0,s)
  
  if(n0 == 0 || n1 == 0) 
  {
    ini <- c(gama_ini,beta_ini,sigma_ini)
  }else{
    pi_ini = rep(0,u)
    ini <- c(gama_ini, pi_ini,beta_ini,sigma_ini)
  }

  #############
  
  opt <- optim(ini, loglik,  escore, # hessian = T, 
                method = "BFGS", control = list(fnscale = -1)) # , maxit = 500, reltol = 1e-9))
  
  if (opt$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")
  
  coef <- opt$par
  k$coef <- coef
  k$conv <- opt$conv
  k$loglik <- opt$value
  k$counts <- as.numeric(opt$counts[1])
  
  k$gamma <- as.vector(coef[1:m])
  k$beta <- as.vector(coef[(m+u+1):(m+u+r)])
  k$sigma <- as.vector(coef[(m+u+r+1):(m+u+r+s)])
  
  eta_hat1 <- as.vector(z1%*%k$gamma)
  eta_hat3 <- as.vector(x1%*%k$beta)
  eta_hat4 <- as.vector(q1%*%k$sigma)
  
  k$eta1 <- eta_hat1
  k$eta3 <- eta_hat3
  k$eta4 <- eta_hat4
  
  k$lambda <- as.vector(linkinv1(eta_hat1))
  k$omega <- as.vector(linkinv3(eta_hat3))
  k$phi <- as.vector(linkinv4(eta_hat4))
  
  if(u!=0)
  {
    k$pi <- as.vector(coef[(m+1):(m+u)])
    k$eta2 <- as.vector(w1%*%k$pi)
    k$p <- as.vector(linkinv2(k$eta2))
  }
  
  # Expected information matrix

  T1 = diag(as.vector( 1/diflink1(k$lambda) ))
  T3 = diag(as.vector( 1/diflink3(k$omega) ))
  T4 = diag(as.vector( 1/diflink4(k$phi) ))

  h1 <- k$omega^(k$phi-2)/ ((1-k$omega^k$phi)*log(1-k$omega^k$phi))
  h2 <- k$omega^(2*k$phi-2)/ ((1-k$omega^k$phi)^(2)*log(1-k$omega^k$phi)^(2))
  delta <- log(0.5)/log(1-k$omega^k$phi)

  kapa<- 0.5772156649
  k0<- (pi^2)/6+kapa^2-2*kapa

  ci = ((k$omega^(k$phi-1))/((1-k$omega^k$phi)*log(1-k$omega^k$phi)))*(delta*log(1-y^k$phi)+1)
  
    b1 <- delta*k$omega^2*h2*log(k$omega)^2*log(1-y^k$phi)
    b2 <- 2*delta*k$omega^2*log(k$omega)*h1*((y^k$phi * log(y))/(1-y^k$phi))
    b3 <- (delta -1)*((y^k$phi * log(y)^2)/((1-y^k$phi)^2))
    b4 <- (ci*k$omega*log(k$omega)^2)*(h1*k$omega^2+ 1/(1-k$omega^k$phi))

    b1 <- as.vector(-1/(k$phi^2) + b1 - b2 - b3 + b4)

    si <- ifelse((i01 == 1),0,b1)

  mi = (k$lambda-1)*k$phi*k$omega*(log(k$omega)*h2+delta*h1*((1-digamma(delta+1)-kapa)/((delta-1)*k$phi)) )

  L = diag(-1/(k$lambda*(1-k$lambda)))
  V = diag((k$lambda-1)*k$phi^2*h2)
  M = diag(mi)
  S = diag(si)
  
  Igg = -t(z1) %*% L %*% (T1^2) %*% z1
  Ibb = -t(x1) %*% V %*% (T3^2) %*% x1
  Ibs = -t(x1) %*% T3 %*% M %*% T4 %*% q1
  Isb = t(Ibs)
  Iss = -t(q1) %*% S %*% (T4^2) %*% q1
  
  k$p <-p
  
  if(u!=0)
  {
    k$pi <- as.vector(coef[(m+1):(m+u)])
    eta_hat2 <- as.vector(w1%*%k$pi)
    k$eta2 <- eta_hat2
    k$p <- as.vector(linkinv2(eta_hat2))
    T2 = diag(as.vector( 1/diflink2(k$p) ))
    P = diag(-k$lambda/(k$p*(1-k$p)))
    Ipp = -t(w1) %*% P %*% (T2^2) %*% w1
    
    I <- rbind(
      cbind(Igg,matrix(0,nrow=m,ncol=u),matrix(0,nrow=m,ncol=r),matrix(0,nrow=m,ncol=s)),
      cbind(matrix(0,nrow=u,ncol=m),Ipp,matrix(0,nrow=u,ncol=r),matrix(0,nrow=u,ncol=s)),
      cbind(matrix(0,nrow=r,ncol=m),matrix(0,nrow=r,ncol=u),Ibb,Ibs),
      cbind(matrix(0,nrow=s,ncol=m),matrix(0,nrow=s,ncol=u),Isb,Iss)
    )
    
  }else{
    I <- rbind(
      cbind(Igg,matrix(0,nrow=m,ncol=r),matrix(0,nrow=m,ncol=s)),
      cbind(matrix(0,nrow=r,ncol=m),Ibb,Ibs),
      cbind(matrix(0,nrow=s,ncol=m),Isb,Iss)
    )
  }
  
  k$I <- I

  ###########################

  vcov <- try(solve(k$I))
  
  k$cov_ok = 0
  if (class(vcov) == "try-error") 
  {
    library(rootSolve)

    k$Knum<-gradient(escore, coef)
    
    vcov <- try(solve(-k$Knum))  # it uses numeric hessian if we have problems with inversion 
    k$cov_ok = 1
    
    if (class(vcov) == "try-error") # 
    {
      vcov <- k$Knum^2 
      warning("Problem with the variance and covariance matrix. The p-values are not correct.")
      k$cov_ok = 2
    }
  }
  
  k$vcov <- vcov
  
  
  loglik_null <- function(theta) 
  {
    lambda <- theta[1]
    p <- theta[2]
    omega <- theta[3]
    phi <- theta[4]
    
    l2 = 0
    if(u!=0)
    {
      l2a = ifelse((y == 1), log(p), 0)
      l2b = ifelse((y == 0), log(1-p), 0)
      l2 = l2a + l2b
    }
    
    ########
    
    l1 = ifelse((i01 == 1), log(lambda), log(1-lambda))
    
    l3 = ifelse((i01 == 1), 0, 
                log(phi)+
                  log(log(0.5)/log(1-omega^phi))+
                  (phi-1)*log(y)+
                  ( (log(0.5)/log(1-omega^phi))-1)*log(1-y^phi)
    )
    
    sum(l1 + l2 + l3)
  }
  
  ini_null<- c(mean(k$lambda),mean(k$p),mean(k$omega),mean(k$phi))
    
  opt_null <- optim(ini_null, loglik_null,  #escore, # hessian = T, 
                    method = "BFGS", control = list(fnscale = -1)) # , maxit = 500, reltol = 1e-9))
  
  
  #####################################
  ### RESIDUALS
  
  # pearson ordinario
  resid1 = (y - k$omega)
  
  # deviance residuals
  l2 = 0
  if(u!=0)
  {
    l2a = ifelse((y == 1), log(k$p), 0)
    l2b = ifelse((y == 0), log(1-k$p), 0)
    l2 = l2a + l2b
  }
  
  l1 = ifelse((i01 == 1), log(k$lambda ), log(1-k$lambda ))
  
  l3 = ifelse((i01 == 1), 0, 
              log(k$phi)+
                log(log(0.5)/log(1-k$omega^k$phi))+
                (k$phi-1)*log(y)+
                ( (log(0.5)/log(1-k$omega^k$phi))-1)*log(1-y^k$phi)
  )
  
  loglikomega <- l1 + l2 + l3
  
  l2 = 0
  if(u!=0)
  {
    l2a = ifelse((y == 1), log(k$p), 0)
    l2b = ifelse((y == 0), log(1-k$p), 0)
    l2 = l2a + l2b
  }
  
  l1 = ifelse((i01 == 1), log(k$lambda ), log(1-k$lambda ))
  
  l3 = ifelse((i01 == 1), 0, 
              log(k$phi)+
                log(log(0.5)/log(1-y^k$phi))+
                (k$phi-1)*log(y)+
                ( (log(0.5)/log(1-y^k$phi))-1)*log(1-y^k$phi)
  )
  
  logliky <- l1 + l2 + l3
  
  resid2 <- sign(y-k$omega)*sqrt(2*(logliky-loglikomega) ) 
  
  # randomized quantile  residuals
  ai0 <- 0 
  bi0 <- k$lambda*(1-k$p)
  
  ai1 <- 1 - k$lambda*k$p
  bi1 <- 1
  
  ru0 <- runif(n,ai0,bi0)
  ru1 <- runif(n,ai1,bi1)

  ui <- (1-k$lambda)*pikum(y,k$lambda,k$p,k$omega,k$phi)
  ui <- ifelse(y==0,ru0,ui)
  ui <- ifelse(y==1,ru1,ui)
  
  resid3 = qnorm(ui)
  
  
  # randomized quantile residuals (other implementation)
  ui<-rep(NA,n)
  for(ci in 1:n)
  {
    if(y[ci]==0) ui[ci] <- runif(1,0,k$lambda[ci]*(1-k$p[ci]))
    if(y[ci]==1) ui[ci] <- runif(1,1-k$lambda[ci],1)
    if(y[ci]!=0 && y[ci]!=1) ui[ci] <- (1-k$lambda[ci])*pikum(y[ci],k$lambda[ci],k$p[ci],k$omega[ci],k$phi[ci])
  }
  resid4 = qnorm(ui)
  
  ######
  k$resid1 <- resid1
  k$resid2 <- resid2
  k$resid3 <- resid3
  k$resid4 <- resid4
  
  #####################################
  
  stderror <- sqrt(diag(k$vcov))
  k$stderror <- stderror
  
  k$zstat <- ((k$coef-theta0)/stderror)
  k$pvalues <- 2*(1 - pnorm(abs(k$zstat) ))
  
  R2 <- 1-exp((-2/n)*(opt$value-opt_null$value))
  k$R2 <- R2
  
  k$aic <- -2* k$loglik + 2*(m+r+s)
  k$bic <- -2* k$loglik+log(n)*(m+r+s)
  
  model_presentation <- cbind(round(k$coef,4),round(k$stderror,4),round(k$zstat,4),round(k$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  k$model <- model_presentation
  
  if(diag==1)
  {
    print(k$model)
  }
  
  return(k)
}