library(tidyverse) 
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")
anole2 <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  print()

#1
anole.log <- anole2%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)
anole.log

#2
anole.lm1 <- lm(HTotal~SVL+ArbPD,data=anole.log)
anole.lm2 <- lm(HTotal~SVL+PH,data=anole.log)


#3:
mod.tibble <- anole.log %>%
  mutate(PH.res = residuals(lmPH)) %>%
  mutate(PD.res = residuals(lmPD)) %>%
  print() 

PH.res.plot <- mod.tibble %>%
  ggplot(aes(PH,PH.res)) + geom_point() + 
  geom_abline(slope=0, col="red") + ylab("Residuals") + xlab("Perch Height") +
  ggtitle("Residuals from Linear Model of Hindlimb-SVL 
          Relationship with Perch Height as Covariate")
PH.res.plot

PD.res.plot <- mod.tibble %>%
  ggplot(aes(ArbPD,PD.res)) + geom_point() + 
  geom_abline(slope=0, col="red") + ylab("Residuals") + xlab("Perch Diameter") +
  ggtitle("Residuals from Linear Model of Hindlimb-SVL 
          Relationship with Perch Diameter as Covariate")
PD.res.plot


#4: 
anole.tree <- read.tree("anole.tre") #to use for phylogenetically correction
pgls.BM.PH <- gls(SVL~HTotal+PH, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML") 
pgls.BM.PD <- gls(SVL~HTotal+ArbPD, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML")
pgls.BM.PH_PD <- gls(SVL~HTotal+PH+ArbPD, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML")


#5: 
anole.phylo.aic <- AICc(pgls.BM.PH, pgls.BM.PD, pgls.BM.PH_PD)
aicw(anole.phylo.aic$AICc)  

anova(pgls.BM.PH) #p-value for PH is 0.0081 (significant for alpha of 0.5 or 0.01): suggests significant predictor of hindlimb length
anova(pgls.BM.PD) #p-value for ArbPD is >0.05 - not a significant predictor
anova(pgls.BM.PH_PD) #p-values for PH and ArbPD are both <0.05: suggests significant predictors

