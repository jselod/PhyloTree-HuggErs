library(tidyverse)
library(ape)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)
library(ggplot2)

#1 establish the data tibble
anole <- read_csv("anole.dat.csv")
anole.eco <- read_csv("anole.eco.csv")
anole.log <- anole%>%
  left_join(anole.eco)%>%
  filter(!Ecomorph%in%c("U","CH"))%>%
  na.omit()%>%
  mutate_at(c("SVL", "HTotal","PH","ArbPD"),log)%>%
  print()


#2 construct two linear models
lm.PH <- lm(SVL~PH+ArbPD,anole.log)
lm.PH
lm.PD <- lm(SVL~ArbPD+PH,anole.log)
lm.PD

#3 residual plots against linear models
#Calculate residuals
anole.log <- anole.log %>%
  mutate(res_PD=residuals(lm.PD))
anole.log <- anole.log %>%
  mutate(res_PH=residuals(lm.PH))
View(anole.log)
#plot residuals
anole.log%>%
  ggplot(aes(ArbPD,res_PD))+
  geom_point() +
  labs(x = "Perch Diameter (PD)", y = "Residuals (PD model)") +
  ggtitle("Residual Plot for PD model")
anole.log%>%
  ggplot(aes(PH,res_PH))+
  geom_point() +
  labs(x = "Perch Height (PH)", y = "Residuals (PH model)") +
  ggtitle("Residual Plot for PH model")


#Q4 3 unique phylogenetic least squares models of the hindlimb-SVL relationships
#PGLS under BM + perch height
anole.tree <- read.tree("anole.tre")
pgls.BM.PH <- gls(HTotal ~SVL + PH, corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")
#PGLS under BM + perch diameter
pgls.BM.PD <- gls(HTotal ~SVL + ArbPD, corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")
#PGLS under BM + perch diameter + perch height + perch diameter
pgls.BM.both <- gls(HTotal ~SVL + PH + ArbPD, corBrownian(1,phy = anole.tree,form=~Species),data = anole.log, method = "ML")

#5 Assess the fit of each of these three models using AICc and AICw 
#AICc from the MuMIn package
anole.phylo.aic <- AICc(pgls.BM.PH,pgls.BM.PD,pgls.BM.both)
aicw(anole.phylo.aic$AICc)
#From the AIC results, a phylogenetically corrected regression model that includes both covariates with traits evolving under BM is the best fit. 
anova(pgls.BM.PH)
anova(pgls.BM.PD)
anova(pgls.BM.both)
#using Anova table, model using perch diameter and model using both covariates are a significant predictor of the hingimb length in a phylogenetic context.

#6 Plot best fitting PGLS model (covariates vs. factors)
anole.log <- anole.log%>%
  mutate(both.res=residuals(pgls.BM.both))
p.eco.phylo <- anole.log%>%
  ggplot(aes(x=Ecomorph,y=both.res)) +geom_boxplot() +stat_summary(fun=mean, geom="point", size=3)
print(p.eco.phylo)



