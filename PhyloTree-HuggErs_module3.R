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
anole.lm1 <- lm(HTotal~ArbPD+Species,data=anole.log)
anole.lm2 <- lm(HTotal~PH+Species,data=anole.log)

#3
anole.resid.1 <- resid(anole.lm1)
anole.resid.2 <- resid(anole.lm2)
