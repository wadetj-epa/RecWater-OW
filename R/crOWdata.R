#create dataset for OW to provide ICF
rm(list=ls())

library(dplyr)
library(timsRstuff)
library(haven)


health<-read_dta("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.dta")
gmeans<-read_dta("L:/Lab/Nheerl_HSF_Beaches/Tim/Rec_Water/EPA_Studies/data/PublicDatasets/geomeans_entero.dta")

#export as text file

write.csv(health, "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.csv", row.names=FALSE)
write.csv(gmeans, "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/geomeans_entero.csv", row.names=FALSE)

health<-read.csv("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/healthfinal.csv", header=TRUE)
gmeans<-read.csv("C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/geomeans_entero.csv",  header=TRUE)
gmeans$date<-as.Date(gmeans$date)
gmeans$strdate<-NULL

#ancillary data
ancil_fresh<-read.csv("L:/Lab/NHEERL_TWade/Tim/Rec_water/EPA_studies/data/PublicDatasets/Archive/Datafiles/Fresh/Ancillary/ancilfresh.txt", header=TRUE)
ancil_marine<-read.csv("L:/Lab/NHEERL_TWade/Tim/Rec_water/EPA_studies/data/PublicDatasets/Archive/Datafiles/Marine/Ancillary/ancilmarine.txt", header=TRUE)
ancil_marine[, c(82, 83)]<-NULL
ancil<-rbind.data.frame(ancil_fresh, ancil_marine)
ancil2<-ancil%>%
  select(c("collectiondate", "beach", "meanairtemp", "meanbirds", "meandogs", 
           "meanbathers",  "cloud8",  "debris8",  "meanbathers",  "rain8",  
           "algae8",  "winddirection8",  "meanwindspeed",  "meanwave", "tide8", 
            "tide15",  "meanuv",  "meanwatertemp" ))

ancil2$date<-as.Date(ancil2$collectiondate, "%m/%d/%Y")

#keep selected health file variables 
health2<-health %>%
  select(c(indid, intdate, sex,  hh ,hcgi, diarrhea, gibase, vomitbase, venfest, anycontact, bodycontact, beach,  age, sex , racecat, swim1face, animunk_any,  gicontact_any, rawmeat_any, comecat, chronany, milescat, bsand, dig, headunder, swallwater))

health2$date<-as.Date(health2$intdate)

drops<-grep("^BB$|^MB$", health2$beach, invert=TRUE)
health2<-health2[drops, ]

health_gmean<-merge(health2, gmeans, by=c("beach", "date"), all=TRUE)
health_gmean<-health_gmean[!is.na(health_gmean$indid), ]

health_gmean_ancil<-merge(health_gmean, ancil2, by=c("beach", "date"), all=TRUE)
health_gmean_ancil<-health_gmean_ancil[!is.na(health_gmean_ancil$indid), ]

healthall<-health_gmean_ancil
#save missing as "." instead of NA
write.csv(healthall, "C:/Users/twade/OneDrive - Environmental Protection Agency (EPA)/Rec_Water/neearcoreall.csv", na=".", row.names=FALSE)


health_mod<-healthall %>%
  filter(venfest!=1 | is.na(venfest)) %>%
  filter(gibase==0|is.na(gibase)) %>%
  filter(vomitbase==0|is.na(vomitbase))


health_mod$exp<-health_mod$log10mean_epcr
health_mod$exp<-ifelse(health_mod$bodycontact==0, 0, health_mod$exp)
health_mod$swimt<-health_mod$bodycontact

otemp<-glm(hcgi~exp+swimt+factor(beach)+age+milescat, data=health_mod, family=binomial(link="identity"))
summary(otemp)
glmCIs(otemp)

#under 10

otemp<-glm(hcgi~exp+swimt+factor(beach)+milescat+meanbathers, data=health_mod, family=binomial(link="identity"), subset=(age<=10))
summary(otemp)
glmCIs(otemp)

#over 10

otemp<-glm(hcgi~exp+swimt+factor(beach)+chronany+milescat+comecat+sex+rawmeat_any+animunk_any, data=health_mod, family=binomial(link="identity"), subset=(age>10))
summary(otemp)
glmCIs(otemp)
