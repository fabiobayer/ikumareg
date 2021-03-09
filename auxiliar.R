# Author: Fabio M Bayer
# Data: March 08, 2021
# e-mail: bayer@ufsm.br


linkfun1 <- function(mu) {
  ret <- log(mu/(1-mu))
  as.vector(ret)
}

linkinv1 <- function(eta) {
  ret <- exp(eta)/(1+exp(eta))
  as.vector(ret)
} 

diflink1 <- function(mu) {
  ret <- 1/(mu*(1-mu))
  as.vector(ret)
}

###### 
linkfun2 <- function(mu) {
  ret <- log(mu/(1-mu))
  as.vector(ret)
}

linkinv2 <- function(eta) {
  ret <- exp(eta)/(1+exp(eta))
  as.vector(ret)
} 

diflink2 <- function(mu) {
  ret <- 1/(mu*(1-mu))
  as.vector(ret)
}


###### 
linkfun3 <- function(mu) {
  ret <- log(mu/(1-mu))
  as.vector(ret)
}

linkinv3 <- function(eta) {
  ret <- exp(eta)/(1+exp(eta))
  as.vector(ret)
} 

diflink3 <- function(mu) {
  ret <- 1/(mu*(1-mu))
  as.vector(ret)
}


###### 
linkfun4 <- function(mu) {
  ret <- log(mu)
  as.vector(ret)
}

linkinv4 <- function(eta) {
  ret <- exp(eta)
  as.vector(ret)
} 

diflink4 <- function(mu) {
  ret <- 1/(mu)
  as.vector(ret)
}

envelope <- function(fit,sim=50,alpha=0.1,name.out=0,name.fig=0)
{
  n <- length(fit$y)
  m <- matrix(0,ncol=n,nrow=sim)

  for(i in 1:sim)
  {
    y = rikum(n,lambda = fit$lambda, p=fit$p, omega = fit$omega, phi = fit$phi)
    fit2 <- ikr.fit(fit$z,fit$w,fit$x,fit$q,y,diag=0)
    m[i,] <- sort(fit2$resid4)
  }
  
  li <- apply(m,2,quantile,alpha/2,na.rm=T)
  med <- apply(m,2,quantile,0.5,na.rm=T)
  ls <- apply(m,2,quantile,1-alpha/2,na.rm=T)
  
  faixa <- range(fit$resid4,li,ls)
  par(pty="s")
  qqnorm(fit$resid4,xlab="Normal quantile",
         ylab="Randomized quantile residuals", ylim=faixa, pch=3)
  par(new=T)
  #
  qqnorm(li,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(ls,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
  par(new=T)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2)
  
  forai<-sum(sort(fit$resid4)<li)
  foras<-sum(sort(fit$resid4)>ls)
  
  fora<-forai+foras
  
  print(c(forai,foras,fora))
  
  
  r<-fit$resid4
  if(name.out!=0) write.table(cbind(r,li,med,ls),file=name.out,quote = F,row.names =F)
  if(name.fig!=0)
  {
    name_pdf<-paste(name.fig,".pdf",sep="")
    pdf(file = name_pdf,width = 4.5, height = 4.2,family = "Times")
    {
      par(mar=c(2.7, 2.7, 0.5, 0.5)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.7, 0))
      
      qqnorm(fit$resid4,xlab="Normal quantile",
             ylab="Randomized quantile residuals", ylim=faixa, pch=3,main=" ")
      par(new=T)
      #
      qqnorm(li,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1,main=" ")
      par(new=T)
      qqnorm(ls,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1,main=" ")
      par(new=T)
      qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2,main=" ")
    }
    dev.off()
    
    name_ps<-paste(name.fig,".eps",sep="")
    postscript(file = name_ps,width = 4.5, height = 4.2,family = "Times",horizontal = FALSE, onefile = FALSE, paper = "special")
    {
      par(mar=c(2.7, 2.7, 0.5, 0.5)) # margens c(baixo,esq,cima,direia)
      par(mgp=c(1.7, 0.7, 0))
      
      qqnorm(fit$resid4,xlab="Normal quantile",
             ylab="Randomized quantile residuals", ylim=faixa, pch=3,main=" ")
      par(new=T)
      #
      qqnorm(li,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1,main=" ")
      par(new=T)
      qqnorm(ls,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1,main=" ")
      par(new=T)
      qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2,main=" ")
    }
    dev.off()
  }  
    
    
}
