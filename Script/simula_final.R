
#Installing and loading Packages
#Install packages

#Packages to be used
packages<-c("readxl","here","tidyverse","ggplot2","gridExtra","knitr","BRugs","coda","rjags","rgl")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Mean function
mea<-function(m,s){exp(m+s^2/2)/(1+exp(m+s^2/2))}

##########################
# Modeling knife to pork #
##########################


#Loading Knife to pork dataset
rep1.1<-c(13,18,12,12,16,38,16,12,9,52,23,14,14,6,12,4,3,12,27,14,21,16,11,3,4,2,3,4,3,2,2,2,3,4,4,3,2,4,2,1,5,10,4,2,2,1,0,2,0,2,1,2,0,1,2,1,0,0,1,0,1,0,0,1,4,1,0,2,1)
rep2.1<-c(20,16,14,20,17,20,18,8,12,46,22,10,16,4,12,4,2,8,36,16,21,12,13,3,3,2,4,5,2,3,2,3,4,2,2,1,2,3,2,0,3,7,4,2,2,1,2,1,0,0,0,1,2,0,0,0,0,0,1,0,0,0,0,0,4,0,1,0,0)
rep3.1<-c(18,16,16,16,16,16,18,8,8,41,24,12,17,6,14,2,2,7,34,16,18,12,15,3,3,3,3,2,2,3,2,2,4,3,2,2,2,3,1,0,2,11,2,2,2,3,1,1,1,1,0,1,2,0,0,0,0,0,0,0,0,0,0,0,9,0,0,1,0)
Y1 = array(cbind(rep1.1, rep2.1, rep3.1),dim=c(23,3,3))

#Model

cat(
  'model{
  
   for(i in 1:23){
    for(j in 1:3){
      for(k in 1:3){
        Y[i,j,k] ~ dpois(D[i,j])
      }
      D[i,j]~dgamma(A[i],l[i,j])
      l[i,j]<-(100/d[j])+0.05
    }
    A[i]<-5+N[i]
    N[i]~dbin(TR[i],tot)
    
    logit(TR[i]) <- beta0[i] 
    beta0[i] ~ dnorm(mu, tau)
  
  }
  tot~dpois(568492)
  
  mu~dnorm(0, 1.0E-6)
  
  tau~dgamma(0.001,0.001)
  
  sigma <-1/(tau)
  
  
  }' , file={f1<-tempfile()} )



# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 10000;nsim = 2000000;nthin = 100

# Defining the data as a list
dadosjags =  list(d=c(10^(-1), 10^(-2),10^(-3)), Y=Y1 )

parms = c("mu", "sigma")

# initializing for adoptation
m1 <- jags.model(f1, dadosjags,  n.chains=3,n.adapt=1000)


# updating burn-in
update(m1,burn)

# Final sampling
mcmc1 <- coda.samples(m1, parms, n.iter=nsim,thin=nthin)

plot(mcmc1)
effectiveSize(mcmc1)
autocorr.diag(mcmc1)
summary(mcmc1)
gelman.plot(mcmc1)
lis1<-summary(mcmc1)

## Result

summary(1/(1+exp(-mcmc1[[1]][,1])))

kable(cbind(round(lis1$statistics[,1:2],2),round(lis1$quantiles,2)))

mea(-3.6,0.4)


linhas<-list()
linhas[[1]]<-(rbinom(100000,5*10^5,  mu1)) 

for (i in 2:50){
linhas[[i]]<-(rbinom(100000,5*10^5,  1/(1+exp(-rnorm(1,lis$statistics[1],lis$statistics[2])))  ))
}

lim1<-(mu1-((1/(1+exp(-lis1$quantiles[1,1])))))*5*10^5
lim2<-((1/(1+exp(-lis1$quantiles[1,5])))-mu1)*5*10^5

plot(density(linhas[[1]]),xlim=c(min(linhas[[1]])-lim2,max(linhas[[1]])+lim1),ylim=c(0,0.005),
     lty=2,lwd = 4,main="",xlab="Number of transfered cells",ylab="",yaxt='n', ann=T,cex=5)

for (i in 2:50){
  if( max(linhas[[i]]) < max(linhas[[1]])+lim2  &  min(linhas[[i]]) > min(linhas[[1]])-lim2  ) {
  
  lines((density(linhas[[i]])))}else {
    NULL
  }

  
}



##########################
# Modeling pork to knife #
##########################


#Loading pork to knife dataset

rep1.2<-c(78,78,95,118,33,75,35,123,122,115,118,118,112,127,27,36,36,47,48,46,51,52,19,19,10,15,3,6,7,8,10,12,16,7,17,16,5,6,4,4,4,3,3,3,1,1,1,3,0,1,1,5,3,0,4,2,3,3,0,2,1,2,2,2,0,1)
rep2.2<-c(88,88,101,142,33,90,59,122,123,120,123,126,137,143,22,41,34,41,45,43,51,50,4,8,15,9,4,5,5,9,9,9,12,12,12,12,4,4,6,3,2,2,4,3,2,2,3,2,0,0,2,0,1,4,2,3,2,2,0,2,1,1,1,1,1,2)
rep3.2<-c(106,106,120,109,22,70,48,58,118,123,111,118,126,137,30,38,33,42,43,43,53,53,8,6,19,19,1,11,4,16,10,6,16,8,15,15,3,4,3,5,2,2,3,2,2,2,2,1,0,1,0,0,2,3,1,2,2,2,0,0,0,1,0,1,0,0)
Y2 = array(cbind(rep1.2, rep2.2, rep3.2),dim=c(22,3,3))

#Model
cat(
  'model{
  
    for(i in 1:22){
    for(j in 1:3){
      for(k in 1:3){
        Y[i,j,k] ~ dpois(D[i,j])
      }
      D[i,j]~dgamma(A[i],l[i,j])
      l[i,j]<-(50/d[j])+0.05
    }
    A[i]<-5+N[i]
    N[i]~dbin(TR[i],tot)
    
    logit(TR[i]) <- beta0[i] 
    beta0[i] ~ dnorm(mu, tau)
  }
  tot~dpois(9411765)
  
  mu~dnorm(0, 1.0E-6)
  
  tau~dgamma(0.001,0.001)
  
  sigma <-1/(tau)
  
  
  }' , file={f2<-tempfile()} )


# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 10000;nsim = 2000000;nthin = 200

# Defining the data as a list
dadosjags2 =  list(d=c(1*10^(-1), 1*10^(-2), 1*10^(-3)), Y=Y2)

parms = c("mu", "sigma")

# initializing for adoptation
m2 <- jags.model(f2, dadosjags2,  n.chains=3,n.adapt=1000)


# updating burn-in
update(m2,burn)

# Final sampling
mcmc2 <- coda.samples(m2, parms, n.iter=nsim,thin=nthin)

plot(mcmc2)
effectiveSize(mcmc2)
autocorr.diag(mcmc2)
summary(mcmc2)
gelman.plot(mcmc2)
lis2<-summary(mcmc2)

## Result

mu2<-1/(1+exp(-lis2$statistics[1]))
ci2<-c(1/(1+exp(-lis2$quantiles[1,1])),1/(1+exp(-lis2$quantiles[1,5])))

kable(cbind(round(lis2$statistics[,1:2],2),round(lis2$quantiles,2)))




linhas2<-list()
for (i in 1:500){
  linhas2[[i]]<-(rbinom(100000,9411765,  1/(1+exp(-rnorm(1,lis2$statistics[1],lis2$statistics[2])))  ))
}

lim1<-600
lim2<-900

plot(density(linhas2[[1]]),xlim=c(min(linhas2[[1]])-lim1,max(linhas2[[1]])+lim1),
     lty=2,lwd = 4,main="",xlab="Number of transfered cells",ylab="",yaxt='n', ann=T,cex=5)

for (i in 2:500){
  if( max(linhas2[[i]]) < max(linhas2[[1]])+lim2  &  min(linhas2[[i]]) > min(linhas2[[1]])-lim2  ) {
    
    lines((density(linhas2[[i]])))}else {
      NULL
    }
  
  
}



## Theoretical

# Valor esperado para p
x<-seq(-5,5,len=100)
y<-seq(0,5,len=100)

z <- outer(x,y, function(x,y) exp(x+y^2/2)/(1+exp(x+y^2/2)))
persp(x,y,z)

persp3d(x,y,z,theta = 30, phi = 30, expand = 0.5,xlab = expression(italic(mu)), ylab = expression(italic(sigma)), 
        zlab = expression(paste("E(",italic(theta),")")))

# Para binomial

#x1<-seq(-5,5,len=100)
#y1<-seq(0,5,len=100)

#z1 <- outer(x1,y1, function(x1,y1) sqrt(exp(x1+y1^2/2)/(1+exp(x1+y1^2/2))*(1-(exp(x1+y1^2/2)/(1+exp(x1+y1^2/2))))) )

#persp(x1,y1,z1)

#persp3d(x,y,z1,col="black")


