#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Model selection - Diversity
# Author: Robin Elahi
# Date: 150629
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# get data for analysis
source("./R/bts_dataPrep.R")

library(nlme)
library(arm)
library(AICcmodavg)

# Random effects for LME models
source("./R/bts_lme_models.R")

# make sure the data are ordered by time within time-series (
# necessary for autocorrelation 
head(divDatSub$year0Z)
tail(divDatSub$year0Z)

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
### FULL DATASET
############################################################
############################################################
lmeDat <- divDat

# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))

# relevel prediction
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")
unique(lmeDat$PredictionOld)

# Set up candidate model list
Cand.mod <- list()

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = div ~ year0Z*Scale + 
                       year0Z*PredictionOld +
                       year0Z*Trophic + year0Z*initialDiv + 
                       year0Z*durationZ, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())
summary(Cand.mod[[1]])

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:PredictionOld)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialDiv)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# year only
Cand.mod[[7]] <- lme(fixed = div ~ year0Z, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

# null model
Cand.mod[[8]] <- lme(fixed = div ~ 1, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", 
              "W/O Year:Trophic Level", "W/O Year:Initial H'", 
              "W/O Year:Study Duration", "Year only", "Null model")

#generate AICc table with names
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/divAIC.csv")

#generate AICc table with numbers
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
###including abundance
############################################################
############################################################
lmeDat <- divDatSub

# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))

# relevel prediction
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")
unique(lmeDat$PredictionOld)

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = div ~ year0Z*Scale + 
                       year0Z*PredictionOld +
                       year0Z*Trophic + year0Z*initialDiv + 
                       year0Z*durationZ + abundCS, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())
summary(Cand.mod[[1]])

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:PredictionOld)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialDiv)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# remove abundance
Cand.mod[[7]] <- update(Cand.mod[[1]], ~. - abundCS)

# year and abundance only
Cand.mod[[8]] <- lme(fixed = div ~ year0Z + abundCS, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

# null model
Cand.mod[[9]] <- lme(fixed = div ~ 1, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", 
              "W/O Year:Trophic Level", "W/O Year:Initial H'", "W/O Year:Study Duration",  
              "W/O Abundance", "Year and abundance only", "Null model")

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_text, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/divAIC reducedAbund.csv")

Cand.modAbund <- Cand.mod[1:9]
Cand.modAbund

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
### REDUCED DATASET, WITHOUT ABUNDANCE IN MODELS
############################################################
############################################################
lmeDat <- divDatSub

# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))

# relevel prediction
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")
unique(lmeDat$PredictionOld)

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = div ~ year0Z*Scale + 
                       year0Z*PredictionOld +
                       year0Z*Trophic + year0Z*initialDiv + 
                       year0Z*durationZ, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:PredictionOld)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialDiv)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# year only
Cand.mod[[7]] <- lme(fixed = div ~ year0Z, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

# null model
Cand.mod[[8]] <- lme(fixed = div ~ 1, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1(), weights = varExp())

Cand.mod
Cand.modDivReduced <- Cand.mod

#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", "W/O Year:Trophic Level", 
              "W/O Year:Initial H'", "W/O Year:Study Duration", "Year only", "Null model")

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_text, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/divAICreduced.csv")
