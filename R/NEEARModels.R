## ----setup, include=FALSE------------------------------------------
knitr::opts_chunk$set(echo = TRUE)



## ----prep, message=FALSE, warning=FALSE----------------------------

#load libraries 
library(haven)
library(dplyr)
library(sandwich)
library(lmtest)
library(ggplot2)
library(Statamarkdown)
#read data
neeardat<-read.csv("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/neearcoreall.csv", na.strings=c(".", "", "NA"))




## ----subset--------------------------------------------------------
neear_mod<-neeardat %>%
  filter(venfest!=1 | is.na(venfest)) %>%
  filter(gibase==0|is.na(gibase)) %>%
  filter(vomitbase==0|is.na(vomitbase))



## ----models--------------------------------------------------------
#exp is the Enterococcus qPCR value
#exp is recoded to 0 for non-swimmers and waders w/o body immersion
#swimt is the swimming definition, in this case body immersion

neear_mod$exp<-neear_mod$log10mean_epcr
neear_mod$exp<-ifelse(neear_mod$bodycontact==0, 0, neear_mod$exp)
neear_mod$swimt<-neear_mod$bodycontact
#All subjects
otempall<-glm(hcgi~exp+swimt+factor(beach)+age+milescat, data=neear_mod, family=binomial(link="identity"))
otempall
round(coeftest(otempall), 4)

#10 and under

otempunder10<-glm(hcgi~exp+swimt+factor(beach)+milescat+meanbathers, data=neear_mod, family=binomial(link="identity"), subset=(age<=10))
otempunder10
round(coeftest(otempunder10), 4)

#over 10

otempover10<-glm(hcgi~exp+swimt+factor(beach)+chronany+milescat+comecat+sex+rawmeat_any+animunk_any, data=neear_mod, family=binomial(link="identity"), subset=(age>10))
otempover10
round(coeftest(otempover10), 4)



## ----cluster-------------------------------------------------------
#cluster standard errors and confidence intervals for all subjects model
clust<-vcovCL(otempall, cluster=neear_mod$hh, cadadjust=TRUE)
round(coeftest(otempall, vcov=clust), 4)
round(coefci(otempall, vcov=clust), 4)



## ----ARfunction----------------------------------------------------
#function to estimate ARs
#bswim-coefficient for swimming (B1); bent-coefficient for entero (B2)
#ent=log10 Enterococcus, vswim=variance of B1; vent=variance of B2
#cov=covariance B1 and B2; level=2 sided confidence interval level
 ARs<-function(bswim, bent,  ent, vswim, vent, cov, level=0.95) {
   z<-qnorm(1-((1-level)/2))
   artemp<-bswim+(bent*ent)
   setemp<-sqrt(vswim+((ent**2)*vent)+(2*ent*cov))
   arup<-artemp+z*setemp
   arlo<-artemp-z*setemp
   return(list(ar=artemp, se=setemp, arup=arup, arlo=arlo))
 }




## ----estARs--------------------------------------------------------
#robust cluster variances and covariances
vswim<-clust[3,3]
vent<-clust[2, 2]
cov<-clust[2, 3]
#set entero levels
ents<-c(seq(1, 3.2, 0.05))
#matrix to store data
armat<-matrix(nrow=length(ents), ncol=3, data=NA)

for(i in 1:length(ents)) {
  entt<-ents[i]
  xtemp<-ARs(bswim=otempall$coefficients[3], bent=otempall$coefficients[2], ent=entt, vswim=vswim, vent=vent, cov=cov)
  armat[i, 1]<-xtemp$ar
  armat[i, 2]<-xtemp$arlo
  armat[i, 3]<-xtemp$arup
}

#create dataset and convert to 1000
ardat<-cbind.data.frame(armat, ents)
names(ardat)<-c("ar", "arlo", "arup", "ent")
ardat$entexp<-10^ardat$ent
ardat$ar1000<-ardat$ar*1000
ardat$arlo1000<-ardat$arlo*1000
ardat$arup1000<-ardat$arup*1000




## ----graph---------------------------------------------------------
#set scales
logbreaks<-c(seq(10, 90, 10), seq(100, 1000, 100))
loglab<-c("10", rep("", 8), "100", rep("", 8), "1000")
arplot<-ggplot(ardat, aes(x=entexp, y=ar1000))+geom_line(linewidth=1.05)+
  scale_x_log10(breaks=logbreaks,labels=loglab)+
  scale_y_continuous(breaks=seq(-20, 60, 10))+
  geom_line(aes(x=entexp, y=arup1000), linetype=2, linewidth=1.05)+
  geom_line(aes(x=entexp, y=arlo1000), linetype=2, linewidth=1.05)+
  labs(x="Enterococcus CCE/100ml \n (daily geometric mean)", y="Attributable Risk per 1000")+ theme_bw()
 
arplot



## ----estARs10------------------------------------------------------
#need to subset data and rerun model
neear_mod10<-dplyr::filter(neear_mod, age<=10)
otempunder10<-glm(hcgi~exp+swimt+factor(beach)+milescat+meanbathers, data=neear_mod10, family=binomial(link="identity"))
#cluster variance
clust10<-vcovCL(otempunder10, cluster=neear_mod10$hh, cadadjust=TRUE)

#robust cluster variances and covariances
vswim10<-clust10[3,3]
vent10<-clust10[2, 2]
cov10<-clust10[2, 3]
#set entero levels (same as above)
ents<-c(seq(1, 3.2, 0.05))
#matrix to store data
armat10<-matrix(nrow=length(ents), ncol=3, data=NA)

for(i in 1:length(ents)) {
  entt<-ents[i]
  xtemp<-ARs(bswim=otempunder10$coefficients[3], bent=otempunder10$coefficients[2], ent=entt, vswim=vswim10, vent=vent10, cov=cov10)
  armat10[i, 1]<-xtemp$ar
  armat10[i, 2]<-xtemp$arlo
  armat10[i, 3]<-xtemp$arup
}

#create dataset and convert to 1000
ardat10<-cbind.data.frame(armat10, ents)
names(ardat10)<-c("ar", "arlo", "arup", "ent")
ardat10$entexp<-10^ardat10$ent
ardat10$ar1000<-ardat10$ar*1000
ardat10$arlo1000<-ardat10$arlo*1000
ardat10$arup1000<-ardat10$arup*1000




## ----graph10-------------------------------------------------------

logbreaks<-c(seq(10, 90, 10), seq(100, 1000, 100))
loglab<-c("10", rep("", 8), "100", rep("", 8), "1000")
arplot10<-ggplot(ardat10, aes(x=entexp, y=ar1000))+geom_line(linewidth=1.05)+
  scale_x_log10(breaks=logbreaks,labels=loglab)+
  scale_y_continuous(breaks=seq(-60, 125, 20))+
  geom_line(aes(x=entexp, y=arup1000), linetype=2, linewidth=1.05)+
  geom_line(aes(x=entexp, y=arlo1000), linetype=2, linewidth=1.05)+
  labs(x="Enterococcus CCE/100ml \n (daily geometric mean)", y="Attributable Risk per 1000")+ theme_bw()
 
arplot10



## ----graphcomb-----------------------------------------------------
#combine ar data to make plotting easier
ardat2<-dplyr::select(ardat, c("entexp", "ar1000", "arlo1000", "arup1000"))
ardat102<-dplyr::select(ardat10, c("ar1000", "arlo1000", "arup1000"))
names(ardat102)<-c("ar1000_10", "arup1000_10", "arlo1000_10")
ardatcomb<-cbind.data.frame(ardat2, ardat102)

arcombplot<-ggplot(ardatcomb, aes(x=entexp, y=ar1000))+geom_line(color="darkred", linewidth=1.05)+
  scale_x_log10(breaks=logbreaks,labels=loglab)+
  scale_y_continuous(breaks=seq(-60, 125, 20))+
  geom_line(aes(x=entexp, y=arup1000), linetype=2, color="darkred", linewidth=1.05)+
  geom_line(aes(x=entexp, y=arlo1000), linetype=2, color="darkred", linewidth=1.05)+
  geom_line(aes(x=entexp, y=ar1000_10), color="navyblue", linewidth=1.05)+
  geom_line(aes(x=entexp, y=arup1000_10), linetype=2, color="navyblue", linewidth=1.05)+
  geom_line(aes(x=entexp, y=arlo1000_10), linetype=2, color="navyblue", linewidth=1.05)+
  labs(x="Enterococcus CCE/100ml \n (daily geometric mean)", y="Attributable Risk per 1000")+ theme_bw()
 

arcombplot


