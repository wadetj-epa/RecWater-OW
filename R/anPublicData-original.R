library(haven)
library(dplyr)
library(timsRstuff)
library(sandwich)
library(lmtest)
library(ggplot2)



health<-read_dta("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.dta")
gmeans<-read_dta("L:/Lab/Nheerl_HSF_Beaches/Tim/Rec_Water/EPA_Studies/data/PublicDatasets/geomeans_entero.dta")

#export as text file
write.csv(health, "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.csv", row.names=FALSE)
write.csv(gmeans, "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/geomeans_entero.csv", row.names=FALSE)

health<-read.csv("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.csv", header=TRUE)
gmeans<-read.csv("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/geomeans_entero.csv",  header=TRUE)


#drop 2009 beaches
drops<-grep("^BB$|^MB$", health$beach, invert=TRUE)
health<-health[drops, ]
table(health$beach)

health_gmean<-merge(health, gmeans, by.x=c("beach", "intdate"), by.y=c ("beach", "date"), all=TRUE)
health_gmean<-health_gmean[!is.na(health_gmean$indid), ]

health_gmean$exp<-health_gmean$log10mean_epcr
health_gmean$exp<-ifelse(health_gmean$bodycontact==0, 0, health_gmean$exp)
health_gmean$swimt<-health_gmean$bodycontact

health_mod<-health_gmean %>%
  filter(venfest!=1 | is.na(venfest)) %>%
  filter(gibase==0|is.na(gibase)) %>%
  filter(vomitbase==0|is.na(vomitbase))


otemp<-glm(hcgi~exp+swimt+factor(beach)+age+milescat, data=health_mod, family=binomial(link="identity"))
 summary(otemp)
 
 glmCIs(otemp)
           
 coeftest(otemp, vcov=sandwich)  
 
 vcovCL(otemp, cluster=health_mod$hh)
 #clust<-vcovCL(otemp, cluster=health_mod$hh, type="HC1")
 #clust<-vcovHC(otemp, cluster=health_mod$hh, type="HC0")
 clust<-vcovCL(otemp, cluster=health_mod$hh, cadadjust=TRUE)
 #stata uses finite simpe adjustment- see if it makes estimates equal
 # M/M-1, where M is number of clusters?
 clust_stata<-vcovCL(otemp, cluster=health_mod$hh, cadadjust=TRUE, type="HC0")
 
M<-length(unique(health_mod$hh))
K<-otemp$rank
N<-otemp$df.null+1
#xx<-(M/(M-1))*((N-1)/(N-K))
xx<-(M/(M-1))
clust_stata<-xx*clust_stata

 
 coeftest(otemp, vcov=clust)
 #close to Stata output for clustered standard errors
 coeftest(otemp)
 
 #calculating standard errors and CIs
 
 dim(clust)
 vswim<-clust[3,3]
 vent<-clust[2, 2]
 cov<-clust[2, 3]
 vswim_stata<-clust_stata[3,3]
 vent_stata<-clust_stata[2, 2]
 cov_stata<-clust_stata[2, 3]
 
 SE3<-sqrt(vswim+((3**2)*vent)+(2*3*cov))
 
 AR3<--0.0273+(3*0.0237)
 AR3_up<-AR3+1.96*SE3
 AR3_lo<-AR3-1.96*SE3
 cbind(AR3_lo, AR3, AR3_up)
 SE2<-sqrt(vswim+((2**2)*vent)+(2*2*cov))
 AR2<--0.0273+(2*0.0237)
 AR2_up<-AR2+1.96*SE2
 AR2_lo<-AR2-1.96*SE2
 cbind(AR2_lo, AR2, AR2_up)

 #nearly recreates SEs from original analysis
 ARs<-function(bswim, bent,  ent, vswim, vent, cov, level=0.95) {
   z<-qnorm(1-((1-level)/2))
   artemp<-bswim+(bent*ent)
   setemp<-sqrt(vswim+((ent**2)*vent)+(2*ent*cov))
   arup<-artemp+z*setemp
   arlo<-artemp-z*setemp
   return(list(ar=artemp, se=setemp, arup=arup, arlo=arlo))
 }
 
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2.039666, vswim=vswim, vent=vent, cov=cov)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=3.159161, vswim=vswim, vent=vent, cov=cov)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.154469, vswim=vswim, vent=vent, cov=cov)

 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2.039666, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=3.159161, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.154469, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 
 x1<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2.039666, vswim=vswim, vent=vent, cov=cov)
 x2<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=3.159161, vswim=vswim, vent=vent, cov=cov)
 x3<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.154469, vswim=vswim, vent=vent, cov=cov)
 
 x1s<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2.039666, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 x2s<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=3.159161, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 x3s<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.154469, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 rbind(c(x1$ar, x1$arlo, x1$arup, x1$se),c(x1s$ar, x1s$arlo, x1s$arup, x1s$se), 
       c(x2$ar, x2$arlo, x2$arup, x2$se),c(x2s$ar, x2s$arlo, x2s$arup, x2s$se),
       c(x3$ar, x3$arlo, x3$arup, x3$se),c(x3s$ar, x3s$arlo, x3s$arup, x3s$se))
 
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=2, vswim=vswim, vent=vent, cov=cov)
 
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.5, vswim=vswim_stata, vent=vent_stata, cov=cov_stata)
 ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=1.5, vswim=vswim, vent=vent, cov=cov)
 
 
  
ents<-c(seq(1, 3.2, 0.05))

armat<-matrix(nrow=length(ents), ncol=3, data=NA)

for(i in 1:length(ents)) {
  entt<-ents[i]
  xtemp<-ARs(bswim=otemp$coefficients[3], bent=otemp$coefficients[2], ent=entt, vswim=vswim, vent=vent, cov=cov)
  armat[i, 1]<-xtemp$ar
  armat[i, 2]<-xtemp$arlo
  armat[i, 3]<-xtemp$arup
}

ardat<-cbind.data.frame(armat, ents)
names(ardat)<-c("ar", "arlo", "arup", "ent")
ardat$entexp<-10^ardat$ent
ardat$ar1000<-ardat$ar*1000
ardat$arlo1000<-ardat$arlo*1000
ardat$arup1000<-ardat$arup*1000

logbreaks<-c(seq(10, 90, 10), seq(100, 1000, 100))
loglab<-c("10", rep("", 8), "100", rep("", 8), "1000")
xx<-ggplot(ardat, aes(x=entexp, y=ar1000))+geom_line()+
  scale_x_log10(breaks=logbreaks, labels=loglab)+scale_y_continuous(breaks=seq(-20, 60, 10))+
  geom_line(aes(x=entexp, y=arup1000), linetype=2)+
  geom_line(aes(x=entexp, y=arlo1000), linetype=2)+
  labs(x="Enterococcus CCE/100ml \n (daily geometric mean)", y="Attributable Risk per 1000")
 
xx



## 10 and under
health_mod10<-dplyr::filter(health_mod, age<=10)

otemp10<-glm(hcgi~exp+swimt+factor(beach)+milescat+meanbathers, data=health_mod10, family=binomial(link="identity"))
clust10<-vcovCL(otemp10, cluster=health_mod10$hh, cadadjust=TRUE)
vswim10<-clust10[3,3]
vent10<-clust10[2, 2]
cov10<-clust10[2, 3]
ARs(bswim=otemp10$coefficients[3], bent=otemp10$coefficients[2], ent=2.039666, vswim=vswim10, vent=vent10, cov=cov10)
ARs(bswim=otemp10$coefficients[3], bent=otemp10$coefficients[2], ent=3.159161, vswim=vswim10, vent=vent10, cov=cov10)
ARs(bswim=otemp10$coefficients[3], bent=otemp10$coefficients[2], ent=1.154469, vswim=vswim10, vent=vent10, cov=cov10)

armat10<-matrix(nrow=length(ents), ncol=3, data=NA)

for(i in 1:length(ents)) {
  entt<-ents[i]
  xtemp<-ARs(bswim=otemp10$coefficients[3], bent=otemp10$coefficients[2], ent=entt, vswim=vswim10, vent=vent10, cov=cov10)
  armat10[i, 1]<-xtemp$ar
  armat10[i, 2]<-xtemp$arlo
  armat10[i, 3]<-xtemp$arup
}

ardat10<-cbind.data.frame(armat10, ents)
names(ardat10)<-c("ar", "arlo", "arup", "ent")
ardat10$entexp<-10^ardat10$ent
ardat10$ar1000<-ardat10$ar*1000
ardat10$arlo1000<-ardat10$arlo*1000
ardat10$arup1000<-ardat10$arup*1000

logbreaks<-c(seq(10, 90, 10), seq(100, 1000, 100))
loglab<-c("10", rep("", 8), "100", rep("", 8), "1000")
xx<-ggplot(ardat10, aes(x=entexp, y=ar1000))+geom_line()+
  scale_x_log10(breaks=logbreaks, labels=loglab)+scale_y_continuous(breaks=seq(-60, 125, 20))+
  geom_line(aes(x=entexp, y=arup1000), linetype=2)+
  geom_line(aes(x=entexp, y=arlo1000), linetype=2)+
  labs(x="Enterococcus CCE/100ml \n (daily geometric mean)", y="Attributable Risk per 1000")

xx


  


                
 
 
#hcgi swimtemp exp i.beach age milescat