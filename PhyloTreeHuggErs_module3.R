### Jaan Project Report ###

#assessing if important ecological parameters of perch height and diameter of the perch
  #predict discrete patterns in the hindlimb-SVL relationship

#loading libraries and datasets:
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

anole <- read_csv("/Users/jaans/Desktop/BIOL3140/anole.dat.csv")
anole_eco <- read_csv("/Users/jaans/Desktop/BIOL3140/anole.eco.csv")


#1: Combining code to establish the anole.log data tibble:
anole2 <- anole %>%
  left_join(anole_eco) %>%
  filter(!Ecomorph %in% c("U", "CH")) %>% 
  na.omit() %>%                         
  print() 
anole.log <- anole2 %>%
  mutate_at(c("SVL", "HTotal", "PH", "ArbPD"), log) %>%
  print()


#2: Constructing 2 simple linear models assessing the effect of perch diameter and height:
lmPH <- lm(SVL~HTotal + PH, data=anole.log)
lmPD <- lm(SVL~HTotal + ArbPD, data=anole.log)


#3: Plot the residuals of models against the discrete factors
mod.tibble <- anole.log %>%
  mutate(PH.res = residuals(lmPH)) %>%
  mutate(PD.res = residuals(lmPD)) %>%
  print() #could have also pivoted table to have values and plot in one plot

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


#4: Constructing PGLS models using Brownian motion and varying combinations of covariates
anole.tree <- read.tree("anole.tre") #to use for phylogenetically correction
pgls.BM.PH <- gls(SVL~HTotal+PH, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML") 
pgls.BM.PD <- gls(SVL~HTotal+ArbPD, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML")
pgls.BM.PH_PD <- gls(SVL~HTotal+PH+ArbPD, correlation=corBrownian(1,phy=anole.tree,form=~Species), data=mod.tibble, method="ML")


#5: Assessing the fit of these models:
anole.phylo.aic <- AICc(pgls.BM.PH, pgls.BM.PD, pgls.BM.PH_PD)
aicw(anole.phylo.aic$AICc)  

anova(pgls.BM.PH) #p-value for PH is 0.0081 (significant for alpha of 0.5 or 0.01): suggests significant predictor of hindlimb length
anova(pgls.BM.PD) #p-value for ArbPD is >0.05 - not a significant predictor
anova(pgls.BM.PH_PD) #p-values for PH and ArbPD are both <0.05: suggests significant predictors

##       fit    delta          w
## 1 -79.62044 2.065530 0.25231023
## 2 -75.88600 5.799965 0.03899534
## 3 -81.68597 0.000000 0.70869443
# we can see that the third model, with both PH and ArbPD as covariates, fits the data best


#6: Plotting effects of covariates on the residuals of best-fitting model:
mod.res <- mod.tibble %>%
  mutate(all.res = residuals(pgls.BM.PH_PD)) %>%
  ggplot(aes(Ecomorph2, all.res)) + geom_boxplot() + 
  stat_summary(fun=mean, geom="point", col="red") + ylab("Residuals") + 
  ggtitle("Residuals from Model of Hindlimb-SVL 
          Relationship with Both Perch Height 
          and Perch Diameter as Co-Variates")
mod.res
