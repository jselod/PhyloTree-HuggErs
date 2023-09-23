### Jaan Selod Module 3 Project ###

#central theme: assess the morphological and allometric differences between ecomorphs in Anolis
 #ecomorph = a group of species not sharing a most recent common ancestor that share an ecological
 # niche and show similar behavior in this niche

#loading libraries:
install.packages("ape")
install.packages("nlme")
install.packages("geiger")
install.packages("caper")
install.packages("phytools")
install.packages("viridis")
install.packages("MuMIn")
library(tidyverse)
library(ape)
library(nlme)
library(geiger)
library(caper)
library(phytools)
library(viridis)
library(MuMIn)

#loading data
anole <- read_csv("anole.dat.csv") #file of morphological data from 46 species of anole
 #includes species, mean snout-vent length (SVL = length between snout and vent), mean total
 # hind-limb length (HTotal)
anole_eco <- read_csv("anole.eco.csv")

anole2 <- anole %>%
  left_join(anole_eco) %>%
  filter(!Ecomorph %in% c("U", "CH")) %>% #removing any rows in the tibble for which the Ecomorph 
  na.omit() %>%                           # value was "U" (unique) or "CH" (chamaeleolis) by keeping
  print()                                 # values NOT equal to those /// also ommitted NAs in any columns

#result: added PH, ArbPD, Ecomorph, and Ecomorph2 as columns

anole.log <- anole2 %>%
  mutate_at(c("SVL", "HTotal", "PH", "ArbPD"), log)
#since our data is continuous, the natural log transformation converts them to proportional representations
# the mutate_at() function mutates columns IN PLACE, replacing each that we identify by wrapping in c()
# then the second argument specifies what we do to the columns
# essentially we can now use a ratio scale for the variables

## Does Hind Limb Length Vary With Size? ##

#let's look at the untransformed data:
anole2 %>%
  ggplot(aes(x=SVL, y=HTotal)) + geom_point() + geom_smooth(method = "lm")
# the data looks to have a linear relationship
# appears that HTotal can be modeled as a linear response to SVL:
 # HTotal = aSVL + b (where a is the slope and b the y intercept)

# let's evaluate how good this fit is with establishing a simple linear model using the lm() function:
anole.lm <- lm(HTotal~SVL, anole2) #NB: x~y --> x predicted by y, 2nd argument is where variables are from
coef(anole.lm)
## (Intercept)         SVL 
##    4.005167    0.639719
anole2 %>%
  ggplot(aes(SVL, HTotal)) + 
  geom_point() + geom_abline(slope = coef(anole.lm)[2], 
                             intercept = coef(anole.lm)[1], 
                             col = "blue")

#with out linear model, we can predict HTotal for a range of snout vent lengths (instead of using geom_smooth)
#let's make a tibble with predictions:
SVL2 <- seq(min(anole2$SVL), max(anole2$SVL), 0.1) #ascending order?
pred.lm <- tibble(
  SVL = SVL2,
  H.pred = predict(anole.lm, newdata = data.frame(SVL = SVL2))
)
anole2 %>%
  ggplot(aes(SVL, HTotal)) +
  geom_point() + 
  geom_point(data = pred.lm, aes(SVL, H.pred), col = "blue")
summary(anole.lm)
## Call:
## lm(formula = HTotal ~ SVL, data = anole2)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -11.8127  -4.6841   0.5629   4.0774  12.3424 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  4.00517    1.92790   2.077   0.0436 *  
## SVL          0.63972    0.02585  24.751   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 6.251 on 44 degrees of freedom
## Multiple R-squared:  0.933,  Adjusted R-squared:  0.9315 
## F-statistic: 612.6 on 1 and 44 DF,  p-value: < 2.2e-16

#summary(anole.lm) is useful! 
  #we have estimated response of HTotal to SVT (the slope 0.63972), and it is significant
  #we also have r^2 value = the percentage of the response variable variation that is explained by
  # our linear model. 0.93 --> this model works pretty well

#but what if the data fits another model better?
#we should assess a non-linear allometric model
#remember that an allometric model is a simple case of an exponential equation, where
  # Htotal = a*SVL^b
  # here, a is the INTERCEPT and b is the SCALING EXPONENT

#to express this model, we can use non-linear least squares (nls):
  #least squares analysis is an approach in regression modeling which approximates the relationship
  #by minimizing the sum of squares of the residuals

#we can fit an allometric model with nls by:
anole.allo <- nls(HTotal~a*SVL^b, start=list(b=1, a=1), data=anole2)
summary(anole.allo)
## 
## Formula: HTotal ~ a * SVL^b
## 
## Parameters:
##   Estimate Std. Error t value Pr(>|t|)    
## b  0.91671    0.03246  28.238  < 2e-16 ***
## a  1.00548    0.14919   6.739 2.76e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 6.102 on 44 degrees of freedom
## 
## Number of iterations to convergence: 4 
## Achieved convergence tolerance: 2.871e-06

#note that we needed to specify starting values for parameters in our model with a list specified by start
# we see here that both parameters a and b are significant, so they are likely explanatory parameters

#but which has a better probability of predicting the data we collected?
#in both, the points of our data fit pretty well - but which model fits our data best?
#we don't have nested models, but models from different families - we need to use AIC
#AICc from the MuMIn package
anole.aic <- AICc(anole.lm,anole.allo)
print(anole.aic) 
#aicw from the geiger package
anole.aicw <- aicw(anole.aic$AICc)
print(anole.aicw)
#both indicate that the allometric model is better, but since the difference between the two is less
  #than 4 AIC units, we can say that allometry and isometry are roughly equivalent models


#now let's consider a more complex relationship of the data:
  #let's consider the effect of ecomorph on hindlimb - SVL relationships

#let's use our log-transformed data and visualize hindlimb-SVL relationships for each ecomorph
  #in ggplot:
anole.log %>%
  ggplot(aes(HTotal, SVL, col = Ecomorph2)) + geom_point() + geom_smooth(method = "lm")
#we can see that ecomorph is an important variable explaining the hindlimb-SVL relationship
  #so let's construct a model that includes this relationship:
anole.log.eco.lm <- lm(HTotal~SVL*Ecomorph2, anole.log)
summary(anole.log.eco.lm)
## Call:
## lm(formula = HTotal ~ SVL * Ecomorph2, data = anole.log)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.217726 -0.037365 -0.000631  0.038581  0.166699 
## 
## Coefficients:
##                           Estimate Std. Error t value Pr(>|t|)   
## (Intercept)                 1.4853     1.0542   1.409  0.16794   
## SVL                         0.6236     0.2108   2.959  0.00559 **
## Ecomorph2grass-bush        -3.1628     1.3919  -2.272  0.02951 * 
## Ecomorph2trunk             -4.9438     2.1623  -2.286  0.02859 * 
## Ecomorph2trunk-crown       -2.3199     1.1646  -1.992  0.05445 . 
## Ecomorph2trunk-ground      -3.3900     1.4537  -2.332  0.02576 * 
## Ecomorph2twig              -2.6170     1.2480  -2.097  0.04351 * 
## SVL:Ecomorph2grass-bush     0.7451     0.3229   2.308  0.02724 * 
## SVL:Ecomorph2trunk          1.1876     0.5365   2.214  0.03365 * 
## SVL:Ecomorph2trunk-crown    0.4744     0.2436   1.948  0.05972 . 
## SVL:Ecomorph2trunk-ground   0.7907     0.3239   2.441  0.02001 * 
## SVL:Ecomorph2twig           0.4963     0.2727   1.820  0.07763 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.08786 on 34 degrees of freedom
## Multiple R-squared:  0.9713, Adjusted R-squared:  0.962 
## F-statistic: 104.7 on 11 and 34 DF,  p-value: < 2.2e-1

anova(anole.log.eco.lm)
## Analysis of Variance Table
## 
## Response: HTotal
##               Df Sum Sq Mean Sq   F value    Pr(>F)    
## SVL            1 8.0567  8.0567 1043.7025 < 2.2e-16 ***
## Ecomorph2      5 0.7560  0.1512   19.5880 3.761e-09 ***
## SVL:Ecomorph2  5 0.0772  0.0154    2.0009    0.1035    
## Residuals     34 0.2625  0.0077                        
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#We performed an ANOVA test on our new model --> this is a two-way analysis of covariance!
  #we are assessing the effect of a categorical variable, Ecomorph2, in the context of HTotal covaries
  # with SVL
#the ANOVA tells us we should REJECT the null hypothesis that ecomorph groups do not have separate 
  #hindlimb-SVL relationships (because P = 3.761e-09)
  # This means that ecomorph has a significant effect on the hindlimb-SVL relationship!

#BUT - does adding the Ecomorph2 parameter result in a better fit compared to a model without it?
  #let's test it out:
anole.log.lm <- lm(HTotal~SVL, anole.log)
anova(anole.log.lm)
## Analysis of Variance Table
## 
## Response: HTotal
##           Df Sum Sq Mean Sq F value    Pr(>F)    
## SVL        1 8.0567  8.0567  323.53 < 2.2e-16 ***
## Residuals 44 1.0957  0.0249                      
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anole.log.aic <- AICc(anole.log.lm,anole.log.eco.lm)
aicw(anole.log.aic$AICc)
##         fit    delta            w
## 1 -34.79896 34.93394 2.595324e-08
## 2 -69.73290  0.00000 1.000000e+00

#yes! the model with the added Ecomorph2 parameter results in a much better fit, both in terms of
  #delta AIC and AICw

#now that we have a better model, we can ask just how much each ecomorph varies among anoles in 
  #their hindlimb-SVL relationship
#let's look at residuals:
anole.log <- anole.log %>%
  mutate(res=residuals(anole.log.lm)) #adds a column with residuals
anole.log%>%
  ggplot(aes(Ecomorph2,res))+geom_point()
  #we can see that the model residuals are greatest in the twig and trunk-ground ecomorphs
  p.eco <- anole.log%>%
    ggplot(aes(x=Ecomorph2,y=res)) +geom_boxplot()
  print(p.eco)
  p.eco + geom_boxplot() + stat_summary(fun=mean, geom="point", size=3)
#now we are representing our residuals in boxplot form, with dots for the mean
  

#By now, we've concluded that ecomorph is an important explanatory variable when it comes to anole
  #hindlimb-SVL relationships
  #but we've ignored the fact that they have an evolutionary history!
  #some species are more closely related than others - they are NOT TRULY INDEPENDENT SAMPLES
  
  #are the patterns we've uncovered so far the result of phylogeny, and not ecomorph alone?
  
#Phylogenetic Comparative Methods (PCMs) = any statistical operation that includes the use of phylogenies
  #in regression analysis, phylogenetic least squares regression (PGLS) is the most flexible
  # and widely used PCM procedure
  # In typical least squares regression operations (performed already), the residuals of any 
  # relationship are assumed to be independent and normally distributed. In PGLS, the variance and 
  # covariance of the residuals are structured given a model of trait and the closeness of 
  # relationships in a phylogenetic tree. 

#So to perform PGLS, we establish this covariation matrix with a tree and assume that the 
  # characters have evolved in this tree under some evolutionary model
  anole.tree <- read.tree("anole.tre") #command from ape package, saving the tree as a phylo object
  plot(anole.tree, cex = 0.4)  #cex = character expansion factor (size of text)
  
  #many PCMs assume a Brownian motion model of character evolution = a random-walk or stumble process
  # this resembles a genetic drift model --> BM can be used to model a process of trait
  # evolution under genetic drift
  
  #another way to model character evolution over a tree is the Ornstein-Uhlenbeck (OU) model
  # traits evolving over a tree are 'pulled' toward an optimum --> often assumed to model a 
  # process of stabilizing selection
  
#PGLS under BM, no ecomorph: 
  pgls.BM1 <- gls(HTotal~SVL, correlation=corBrownian(1, phy=anole.tree, form=~Species), data=anole.log, method="ML")
#PGLS under BM, with ecomorph:
  pgls.BM2 <- gls(HTotal~SVL*Ecomorph2, correlation=corBrownian(1, phy=anole.tree, form=~Species), data=anole.log, method="ML")

#PGLS under OU, no ecomorph:
  pgls.OU1 <- gls(HTotal~SVL, correlation=corMartins(0, phy=anole.tree, form=~Species), data=anole.log, method="ML")
#PGLS under OU, with ecomorph:
  pgls.OU2 <- gls(HTotal~SVL*Ecomorph2, correlation=corMartins(0, phy=anole.tree, form=~Species), data=anole.log, method="ML")

#what we did above: set up generalized least squares models that establish a correlation matrix
  #based on BM and OU models of evolution
  #the matrix is established with corBrownian() for BM models and corMartins() for OU models
  # note: we don't need to specify any parameters --> they are estimated under maximum likelihood
  # when method = "ML"
  
#now we can assess which models fit the data best!
  anole.phylo.aic <- AICc(pgls.BM1, pgls.BM2, pgls.OU1, pgls.OU2)
  aicw(anole.phylo.aic$AICc)  
  ##      fit     delta            w
  ## 1 -63.66118 29.319176 4.225281e-07
  ## 2 -92.98036  0.000000 9.827289e-01
  ## 3 -57.31864 35.661714 1.772520e-08
  ## 4 -84.89770  8.082653 1.727062e-02  
  # we can see that the phylogenetically corrected regression model including Ecomorph2 with traits
    # evolving under BM is the best fit
    #we can interpret this to mean that the traits have evolved randomly, but randomly within each lineage
  
#now, we can consider if Ecomorph is a significant factor in predicting the hindlimb-SVL relationship
  #in a phylogenetic context:
  anova(pgls.BM2)
  ## Denom. DF: 34 
  ##               numDF   F-value p-value
  ## (Intercept)       1 15848.498  <.0001
  ## SVL               1  1073.374  <.0001
  ## Ecomorph2         5    16.731  <.0001
  ## SVL:Ecomorph2     5     1.596  0.1878
  # we can see that the p-value for Ecomorph2 is significant for the effect on the relationship
  
#but when we visualize, we can see that ecomorph isn't as strong a factor when we consider phylogeny:
  #let's look at the residuals of the hindlimb-SVL relationship, but now considering BM evolution of 
  # this relationship over the tree
  anole.log <- anole.log %>%
    mutate(phylo.res = residuals(pgls.BM2)) #adding a column for phylogenetically corrected residuals
  p.eco.phylo <- anole.log %>%
    ggplot(aes(x=Ecomorph2, y=phylo.res)) + geom_boxplot() + stat_summary(fun=mean, geom="point", size=3)
  print(p.eco.phylo) #plotting those residuals against ecomorph
  #now that we have the phylogenetically corrected residuals, we can compare them to our previous
  # uncorrected residuals
  #let's make the anole.log tibble longer with respect to the two types of residuals
  anole.log %>%
    dplyr::select(Ecomorph2,res,phylo.res) %>% #selected the 3 columns of interest
    pivot_longer(cols=c("res","phylo.res")) %>% #specified columns to make the tibble longer
    print %>%                                   #now, for each entry of an individual ecomorph, there are 2 types of residuals identified by name and value!
    ggplot(aes(x=Ecomorph2,y=value)) + geom_boxplot() + 
      stat_summary(fun=mean, geom="point", size=3) + facet_grid(name~.,scales="free_y") + ylab("residual")
  #value --> column in our tibble with values from both types of residuals
  #for facet_grid(), name~. referred to how row facets ~ column facets --> we want each of the values
    # according to name (the type of residual)!
  #scales = "free_y" allows the scale on the y axis to be different in each row
  
  #we can now see that the residuals condense when we consider phylogeny - there was a lot of
  # phylogenetic signal in our hindlimb-SVL data
  #when we consider phylogeny, the patterns of which groups vary changes
    #the trunk-ground ecomorph no longer has such a high residual - rather, it is pretty close
    # to the phylogenetically corrected mean for anoles - this is an unremarkable group
  