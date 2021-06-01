#Loading packages and functions

#Packages to be used
packages<-c("here","tidyverse","ggplot2","gridExtra","knitr","lme4")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#Function for back transformation
back<-function(x){1/(1+exp(-x))}



#Loading Knife to pork dataset
rep1.1<-c(13,18,12,12,16,38,16,12,9,52,23,14,14,6,12,4,3,12,27,14,21,16,11,3,4,2,3,4,3,2,2,2,3,4,4,3,2,4,2,1,5,10,4,2,2,1,0,2,0,2,1,2,0,1,2,1,0,0,1,0,1,0,0,1,4,1,0,2,1)
rep2.1<-c(20,16,14,20,17,20,18,8,12,46,22,10,16,4,12,4,2,8,36,16,21,12,13,3,3,2,4,5,2,3,2,3,4,2,2,1,2,3,2,0,3,7,4,2,2,1,2,1,0,0,0,1,2,0,0,0,0,0,1,0,0,0,0,0,4,0,1,0,0)
rep3.1<-c(18,16,16,16,16,16,18,8,8,41,24,12,17,6,14,2,2,7,34,16,18,12,15,3,3,3,3,2,2,3,2,2,4,3,2,2,2,3,1,0,2,11,2,2,2,3,1,1,1,1,0,1,2,0,0,0,0,0,0,0,0,0,0,0,9,0,0,1,0)
Y1 = array(cbind(rep1.1, rep2.1, rep3.1),dim=c(23,3,3))

data1<-cbind.data.frame(
id=c(rep(1:23,9)),

d=c(rep(0.1,23),rep(0.01,23),rep(0.001,23),rep(0.1,23),rep(0.01,23),rep(0.001,23),rep(0.1,23),rep(0.01,23),rep(0.001,23)),

rep=c(rep(1,23*3),rep(2,23*3),rep(3,23*3)),

cfu=c(
rep1.1,
rep2.1,
rep3.1)

)

data1$cfu1<-data1$cfu/data1$d*100

data1$total<-rep(568492,207)
data1$neg<-data1$total-data1$cfu


summary(glmer(cbind(data1$cfu1,data1$neg)~1+(1|d)+(1|id),data=data1,family = binomial(link = "logit"))->total)


#Shrink
c<-15/16*pi/sqrt(3)

sr<-1/sqrt((sum((total@theta)^2)/c^2)+1)



back(summary(total)$coefficients[1]*sr)


n_inc<-1000
inc<-rnorm(n_inc,summary(total)$coefficients[1]*sr,total@theta[2]+total@theta[3])

var<-list()
med1<-list()

for(i in 1:n_inc){
    var[[i]]<-rnorm(10000,inc[i],total@theta[1])
    med1[i]<-mean(var[[i]])
}


plot(ecdf(back(var[[1]])))
for(i in 2:n_inc){
    lines((ecdf(back(var[[i]]))))
  
}

back(mean(unlist(var)))
back(median(unlist(var)))

#Loading pork to knife
rep1.2<-c(78,78,95,118,33,75,35,123,122,115,118,118,112,127,27,36,36,47,48,46,51,52,19,19,10,15,3,6,7,8,10,12,16,7,17,16,5,6,4,4,4,3,3,3,1,1,1,3,0,1,1,5,3,0,4,2,3,3,0,2,1,2,2,2,0,1)
rep2.2<-c(88,88,101,142,33,90,59,122,123,120,123,126,137,143,22,41,34,41,45,43,51,50,4,8,15,9,4,5,5,9,9,9,12,12,12,12,4,4,6,3,2,2,4,3,2,2,3,2,0,0,2,0,1,4,2,3,2,2,0,2,1,1,1,1,1,2)
rep3.2<-c(106,106,120,109,22,70,48,58,118,123,111,118,126,137,30,38,33,42,43,43,53,53,8,6,19,19,1,11,4,16,10,6,16,8,15,15,3,4,3,5,2,2,3,2,2,2,2,1,0,1,0,0,2,3,1,2,2,2,0,0,0,1,0,1,0,0)
Y2 = array(cbind(rep1.2, rep2.2, rep3.2),dim=c(22,3,3))


data2<-cbind.data.frame(
  id=c(rep(1:22,9)),
  
  d=c(rep(0.1,22),rep(0.01,22),rep(0.001,22),rep(0.1,22),rep(0.01,22),rep(0.001,22),rep(0.1,22),rep(0.01,22),rep(0.001,22)),
  
  rep=c(rep(1,22*3),rep(2,22*3),rep(3,22*3)),
  
  cfu=c(
    rep1.2,
    rep2.2,
    rep3.2)
  
)

data2$cfu1<-data2$cfu/data2$d*50

data2$total<-rep(9411765,198)
data2$neg<-data2$total-data2$cfu



summary(glmer(cbind(data2$cfu1,data2$neg)~1+(1|id)+(1|d/rep),data=data2,family = binomial(link = "logit"))->total2)


#Shrink
c<-15/16*pi/sqrt(3)

sr2<-1/sqrt((sum((total2@theta)^2)/c^2)+1)

back(summary(total2)$coefficients[1]*sr2)




inc2<-rnorm(1000,summary(total2)$coefficients[1]*sr,total2@theta[2]+total2@theta[3])

var2<-list()


for(i in 1:1000){
  var2[[i]]<-rnorm(10000,inc2[i],total2@theta[1])
  
}


plot(ecdf(back(var2[[1]]) ))
for(i in 2:100){
  lines((ecdf(back(var2[[i]]))))
}

back(mean(unlist(var2)))
back(median(unlist(var2)))
