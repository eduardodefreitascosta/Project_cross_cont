
#Installing and loading Packages

if(!require(BRugs)){
  install.packages("BRugs")
}
library(BRugs)



if(!require(coda)){
  install.packages("coda")
}
library(coda)


#Loading Knife to pork dataset
rep1.1<-c(13,18,12,12,16,38,16,12,9,52,23,14,14,6,12,4,3,12,27,14,21,16,11,3,4,2,3,4,3,2,2,2,3,4,4,3,2,4,2,1,5,10,4,2,2,1,0,2,0,2,1,2,0,1,2,1,0,0,1,0,1,0,0,1,4,1,0,2,1)
rep2.1<-c(20,16,14,20,17,20,18,8,12,46,22,10,16,4,12,4,2,8,36,16,21,12,13,3,3,2,4,5,2,3,2,3,4,2,2,1,2,3,2,0,3,7,4,2,2,1,2,1,0,0,0,1,2,0,0,0,0,0,1,0,0,0,0,0,4,0,1,0,0)
rep3.1<-c(18,16,16,16,16,16,18,8,8,41,24,12,17,6,14,2,2,7,34,16,18,12,15,3,3,3,3,2,2,3,2,2,4,3,2,2,2,3,1,0,2,11,2,2,2,3,1,1,1,1,0,1,2,0,0,0,0,0,0,0,0,0,0,0,9,0,0,1,0)
Y1 = array(cbind(rep1.1, rep2.1, rep3.1),dim=c(23,3,3))


#Loading pork to knife dataset

rep1.2<-c(78,78,95,118,33,75,35,123,122,115,118,118,112,127,27,36,36,47,48,46,51,52,19,19,10,15,3,6,7,8,10,12,16,7,17,16,5,6,4,4,4,3,3,3,1,1,1,3,0,1,1,5,3,0,4,2,3,3,0,2,1,2,2,2,0,1)
rep2.2<-c(88,88,101,142,33,90,59,122,123,120,123,126,137,143,22,41,34,41,45,43,51,50,4,8,15,9,4,5,5,9,9,9,12,12,12,12,4,4,6,3,2,2,4,3,2,2,3,2,0,0,2,0,1,4,2,3,2,2,0,2,1,1,1,1,1,2)
rep3.2<-c(106,106,120,109,22,70,48,58,118,123,111,118,126,137,30,38,33,42,43,43,53,53,8,6,19,19,1,11,4,16,10,6,16,8,15,15,3,4,3,5,2,2,3,2,2,2,2,1,0,1,0,0,2,3,1,2,2,2,0,0,0,1,0,1,0,0)
Y2 = array(cbind(rep1.2, rep2.2, rep3.2),dim=c(22,3,3))



############################
# Modeling knife para meat #
############################

modelo1 = function(){
  for(i in 1:23){
    for(j in 1:3){
      for(k in 1:3){
        Y[i,j,k] ~ dpois(D[i,j])
      }
      D[i,j]~dgamma(A[i],l[i,j])
      l[i,j]<-(100/d[j])+0.05
    }
    A[i]<-5+N[i]
    N[i]~dbin(TR[i],568492)
    
    logit(TR[i]) <- beta0[i] 
    beta0[i] ~ dnorm(mu, tau)
  }
  
  mu~dunif(-10,10)
  tau <-1/(sigma)
  sigma~dunif(0,1000)
  TR.medio <- exp(mu)/(1+exp(mu))
  ICC <- (sigma*sigma)/(sigma*sigma + 3.289868)
  
}

inits1 <- function(){
  list(sigma = 10, mu=-3.4,beta0=rep( -3.45,23), N=rep(50000,23),D=rep(500,69))
}

parameters1 <- c( "TR.medio", "sigma", "ICC",'mu')   


dados1 = list(d=c(10^(-1), 10^(-2),10^(-3)), Y=Y1 )


simul_faca_carne <- BRugsFit(modelo1, dados1, inits1, numChains = 1, parameters1,
                             nBurnin = 10000, nIter=6000000, nThin=150, coda=T)


plot(simul_faca_carne )
kable(effectiveSize(simul_faca_carne))
kable(autocorr.diag(simul_faca_carne))
summary(simul_faca_carne)


##########################
# Modeling meat to knife #
##########################

modelo2 = function(){
  for(i in 1:22){
    for(j in 1:3){
      for(k in 1:3){
        Y[i,j,k] ~ dpois(D[i,j])
      }
      D[i,j]~dgamma(A[i],l[i,j])
      l[i,j]<-(50/d[j])+0.05
    }
    A[i]<-5+N[i]
    N[i]~dbin(TR[i],9411765)
    
    logit(TR[i]) <- beta0[i] 
    beta0[i] ~ dnorm(mu, tau)
  }
  
  mu~dunif(-10,10)
  tau <-1/(sigma)
  sigma~dunif(0,1000)
  TR.medio <- exp(mu)/(1+exp(mu))
  ICC <- (sigma*sigma)/(sigma*sigma + 3.289868)
  
}



inits2 <- function(){
  list(sigma = 10, mu=-3.4,beta0=rep( -3.45,23), N=rep(50000,23),D=rep(500,69))
}

parameters2 <- c( "beta0")   




dados2 = list(d=c(1*10^(-1), 1*10^(-2), 1*10^(-3)), Y=Y2)

simul_carne_faca <- BRugsFit(modelo2, dados2, inits2, numChains = 1, parameters2,
                             nBurnin = 25000, nIter=6000000, nThin=150, coda=T)


plot(simul_carne_faca)
effectiveSize(simul_carne_faca)
autocorr.diag(simul_carne_faca)
summary(simul_carne_faca)



##Descriptive classical approach 1� dilution

par(mfrow=c(1,2))

freq<-apply(matrix(apply(Y1[1:23,1,],1,mean)/(c(rep(0.1,23))/100),nrow=23,ncol=1),1,mean ) /568492.06

hist(freq, ylab="Frequency",xlab="Transfer ratio",main=" ")
par(xpd=T)
text(0,12.5,labels="A")
mean(freq)
sd(freq)
summary(freq)


freq2<-apply(matrix(apply(Y2[1:22,1,],1,mean)/(c(rep(0.1,22))/50),nrow=22,ncol=1),1,mean ) /9411765

hist(freq2, ylab="Frequency",xlab="Transfer ratio",main=" ")
par(xpd=T)
text(0.001,8,labels="B")
mean(freq2)
sd(freq2)
summary(freq2)



# Calculation of transfer probability
TRs1<- simul_faca_carne[,2,1]
min(TRs1[[1]])
max(TRs1[[1]])
mean(TRs1[[1]])
quantile(TRs1[[1]],0.5)
quantile(TRs1[[1]],0.95)



TRs2<- simul_carne_faca[,2,1]
min(TRs2[[1]])
max(TRs2[[1]])
mean(TRs2[[1]])
quantile(TRs2[[1]],0.5)


# Calculation of the mean sigma
sig1<- simul_faca_carne[,4,1]
sig2<- simul_carne_faca[,4,1]


#Calculation of the mean ICC
ICC1<- simul_faca_carne[,1,1]
ICC2<- simul_carne_faca[,1,1]

#Calculation of the mean mu

mu1<-simul_faca_carne[,3,1]
mu2<-simul_carne_faca[,3,1]



#ECDF and density plots


par(mfrow=c(1,2))
plot(density(TRs1[[1]]),main='',xlab='Transfer probability',lwd = 2,xlim=c(0.01,0.045),
     panel.first = rect(c(0,0.035), 0, c(0.02,10), 1e6, col='grey', border=NA)     )
text(0.01, 121, paste('A'), xpd=NA)
abline(v=quantile(TRs1[[1]],0.025));abline(v=quantile(TRs1[[1]],0.975))


plot(density(TRs2[[1]]),main='',xlab='Transfer probability',lwd = 2,xlim=c(0.002,0.006),
     panel.first = rect(c(0,quantile(TRs2[[1]],0.975)), 0, c(quantile(TRs2[[1]],0.025),10), 1e6, col='grey', border=NA))
text(0.002, 954, paste('B'), xpd=NA)
abline(v=quantile(TRs2[[1]],0.025));abline(v=quantile(TRs2[[1]],0.975))





#Stochastic model with ICC spliting Variability and uncertainty

#Stochastic model knife to meat

nuncer<-10000
nvar<-10000
repl<-100

unc1<-numeric()
simul1<-numeric()
nor1<-numeric()
med_1<-numeric()

par(mfrow=c(1,2))
plot(length(simul1),xlim=c(0,0.2),ylim=c(0,1),main='',xlab='Transfer probability',col='grey')
text(0, 1.1, paste('A'), xpd=NA)


unc1<-rnorm(nuncer, mean(mu1[[1]]),sqrt(mean(sig1[[1]]))*sqrt(1-(mean(ICC1[[1]]))))

for (i in 1:repl){
  nor1<-rnorm(nvar,sample(unc1,size=1,prob=rep((1/nuncer),nuncer),replace=F),sqrt(mean(sig1[[1]]))*sqrt((mean(ICC1[[1]]))))
  simul1<-exp(nor1)/(1+exp(nor1))
  
  
  lines(ecdf(simul1),col='grey')
  med_1[i]<-median(simul1)
}

arrows(x0=quantile(med_1,0.025), y0=0.5, x1=quantile(med_1,0.975), y1=0.5,code = 3,lwd = 2)




#Stochastic model meat ot knife

unc2<-numeric()
simul2<-numeric()
nor2<-numeric()
med_2<-numeric()

plot(length(simul2),xlim=c(0,0.02 ),ylim=c(0,1),main='',xlab='Transfer probability',col='grey')
text(0, 1.1, paste('B'), xpd=NA)




unc2<-rnorm(nuncer, mean(mu2[[1]]),sqrt(mean(sig2[[1]]))*sqrt(1-(mean(ICC2[[1]]))))

for(i in 1:repl){
  nor2<-rnorm(nvar,sample(unc2,size=1,prob=rep((1/nuncer),nuncer),replace=F),sqrt(mean(sig2[[1]]))*sqrt((mean(ICC2[[1]]))))
  simul2<-exp(nor2)/(1+exp(nor2))
  lines(ecdf(simul2),col='grey')
  med_2[i]<-median(simul2)
}


arrows(x0=quantile(med_2,0.025), y0=0.5, x1=quantile(med_2,0.975), y1=0.5,code = 3,lwd = 2)




#Medians intervals
med_1<-apply(simul1,2, median)
med_2<-apply(simul2,2, median)
summary(med_1)
summary(med_2)

media1<-apply(simul1,2, mean)
media2<-apply(simul2,2, mean)


#Min, Max, and percentiles
min(simul1)
max(simul1)
mean(simul1)

min(simul2)
max(simul2)


quantile(simul1,0.025)
quantile(simul1,0.975)

quantile(simul2,0.025)
quantile(simul2,0.975)

#Means
mean(media1)
mean(media2)
