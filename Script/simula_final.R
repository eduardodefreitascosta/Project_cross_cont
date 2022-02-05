
#header #################################################################################
#'simula_final.R'

#Title: Salmonella transfer
#Project ID: pid
#Client: UFRGS
#Author: <Eduardo> <Costa>, Wageningen Bioveterinary Research

#Description: This is a paper describying the a Bayesian model to estimate the *Salmonella* sp. transfer probability between knife and pork in a household scenario.

#Start date: date
#Last Update: {6:date}

#R version: r.version
#Scriptversion: version

#Dependencies
#<-Downstream none
#->Upstream Main.R

#Input:
#- None

#Output:
#- pork_to_knife.txt; knife_to_pork.txt

#Peer reviewer(s)

#Please ensure directories are relative. Please comment code sufficiently.

#Script start#############################################################################



##########################
# Modeling knife to pork #
##########################


#Loading Knife to pork dataset
rep1.1<-c(13,18,12,12,16,38,16,12,9,52,23,14,14,6,12,4,3,12,27,14,21,16,11,3,4,2,3,4,3,2,2,2,3,4,4,3,2,4,2,1,5,10,4,2,2,1,0,2,0,2,1,2,0,1,2,1,0,0,1,0,1,0,0,1,4,1,0,2,1)
rep2.1<-c(20,16,14,20,17,20,18,8,12,46,22,10,16,4,12,4,2,8,36,16,21,12,13,3,3,2,4,5,2,3,2,3,4,2,2,1,2,3,2,0,3,7,4,2,2,1,2,1,0,0,0,1,2,0,0,0,0,0,1,0,0,0,0,0,4,0,1,0,0)
rep3.1<-c(18,16,16,16,16,16,18,8,8,41,24,12,17,6,14,2,2,7,34,16,18,12,15,3,3,3,3,2,2,3,2,2,4,3,2,2,2,3,1,0,2,11,2,2,2,3,1,1,1,1,0,1,2,0,0,0,0,0,0,0,0,0,0,0,9,0,0,1,0)
Y1 = array(cbind(rep1.1, rep2.1, rep3.1),dim=c(23,3,3))

#Knife dataset
rep1.3<-c(82,127,117,133,58,74,62,70,107,129,127,96,144,127,124,126,137,126,128,122,114,6,11,15,13,14,6,12,12,18,10,12,10,16,12,16,13,11,10,5,5,7)
rep1.3<-c(70,164,137,147,70,105,39,62,114,139,122,98,136,129,136,122,139,128,129,120,118,13,18,12,21,5,9,5,10,12,9,19,14,12,16,13,14,11,8,8,5,8)
rep1.3<-c(64,144,124,168,57,121,54,63,100,140,125,124,128,137,141,118,128,131,126,118,121,9,11,17,15,13,13,10,7,10,12,16,12,13,14,15,14,6,7,6,6,12)
Y3 = array(cbind(rep1.3, rep1.3, rep1.3),dim=c(21,2,3))



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
    N[i]~dbin(TR,tot)
  
   }
    TR~dbeta(a,b)
    
    a~dunif(0,1000)
    b~dunif(0,1000)
   
  tot~dpois(568492)
  
  }' , file={f1<-tempfile()} )



# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 100000;nsim = 2000000;nthin = 1000

# Defining the data as a list
dadosjags =  list(d=c(10^(-1), 10^(-2),10^(-3)), Y=Y1,d2=c(10^(-1), 10^(-2)),W=Y3)

parms = c("TR","a","b")

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


capture.output(summary(mcmc1), file = here("Output","knife_to_pork.txt"))

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
    N[i]~dbin(TR,tot)
    
    

  }
    TR ~ dbeta(a, b)  
  tot~dpois(9411765)
  a~dunif(0,1000)
    b~dunif(0,1000)
  
  }' , file={f2<-tempfile()} )




# Defining burn-in, number of simulations, thin and parameters to be shown
burn = 100000;nsim = 2000000;nthin = 1000

# Defining the data as a list
dadosjags2 =  list(d=c(1*10^(-1), 1*10^(-2), 1*10^(-3)), Y=Y2)

parms = c("TR", "a","b")

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

capture.output(summary(mcmc2), file = here("Output","pork_to_knife.txt"))


###############################################
# Descriptive classical approach 1st dilution #
###############################################

par(mfrow=c(1,2))

# Knife to meat
freq<-apply(matrix(apply(Y1[1:23,1,],1,mean)/(c(rep(0.1,23))/100),nrow=23,ncol=1),1,mean ) /568492.06


# Meat to knife
freq2<-apply(matrix(apply(Y2[1:22,1,],1,mean)/(c(rep(0.1,22))/50),nrow=22,ncol=1),1,mean ) /9411765


# Total dataset descriptive
teste<-cbind.data.frame(freq=c(freq,freq2),label=c(rep("K-M",length(freq)),rep("M-K",length(freq2))) )
describeBy(teste$freq,group=teste$label,digits = 5,mat=T)


# Layout to split the screen


tiff(here("Figures","Fig_11.tiff"), units="in", width=6, height=5, res=300)

par(mar=c(0.1, 3, 1.1, 2.1))
layout(mat = (matrix(c(1,2,3,4), byrow=TRUE)),  height = c(1.5,8,1.5,8))

# Draw the boxplot and the histogram 
par(mar=c(0, 2, 1.1, 2.1))
boxplot(freq , horizontal=TRUE , ylim=c(0,max(freq)+0.02), xaxt="n" , frame=F)
par(mar=c(4, 2, 1.1, 2.1))
hist(freq , breaks=10  , border=T , main="", yaxt="n", xlab="Transfer probability", xlim=c(0,max(freq)+0.02),cex.lab=1.5, 
     cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
text(0,10,"A",cex=2)
abline(v=quantile(freq,c(0.025,0.975)))

# Draw the boxplot and the histogram 
par(mar=c(0.1, 2, 1.1, 2.1))
boxplot(freq2 , horizontal=TRUE , ylim=c(0,max(freq2)+0.002), xaxt="n" , frame=F)
par(mar=c(4.5, 2, 1.1, 2.1))
hist(freq2 , breaks=10  , border=T , main="", yaxt="n", xlab="Transfer probability", xlim=c(0,max(freq2)+0.002)
     ,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
text(0,3.8,"B",cex=2)
abline(v=quantile(freq2,c(0.025,0.975)))


dev.off()
