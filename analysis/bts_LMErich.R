#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Model selection - Richness
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

# make sure the data are ordered by time within time-series 
# necessary for autocorrelation 
head(richDat$year0Z)
tail(richDat$year0Z)

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: RICHNESS
### FULL DATASET
############################################################
############################################################
lmeDat <- richDat
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = richLN ~ 1 +  
                       year0Z*Scale + year0Z*Prediction +
                       year0Z*Trophic + year0Z*initialRichLN + 
                       year0Z*durationZ, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1())
summary(Cand.mod[[1]])$tTable

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:Prediction)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialRichLN)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# year only
Cand.mod[[7]] <- lme(fixed = richLN ~ 1 +  
	year0Z, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())

# null model
Cand.mod[[8]] <- lme(fixed = richLN ~ 1, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())

#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", 
              "W/O Year:Trophic Level", "W/O Year:Initial S", "W/O Year:Study Duration", 
              "Year only", "Null model")

#generate AICc table with numbers
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)

#generate AICc table with names
mod.aicctab <- aictab(cand.set= Cand.mod, modnames= mod_text, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
# write.csv(mod.aicctab, "./output/richAIC.csv")

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: RICHNESS
### REDUCED DATASET, WITH ABUNDANCE IN MODELS
############################################################
############################################################

lmeDat <- subDat
unique(lmeDat$Prediction)
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")
dim(lmeDat)
levels(lmeDat$Trophic)

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = richLN ~ 1 +  
	year0Z*Scale + year0Z*Prediction +
	year0Z*Trophic + year0Z*initialRichLN + 
	year0Z*durationZ + abundCS, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())

summary(Cand.mod[[1]])
Cand.mod[[1]]

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:Prediction)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialRichLN)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# remove abundance
Cand.mod[[7]] <- update(Cand.mod[[1]], ~. - abundCS)

# year and abundance only
Cand.mod[[8]] <- lme(fixed = richLN ~ year0Z + abundCS, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1())

# null model
Cand.mod[[9]] <- lme(fixed = richLN ~ 1, 
                     data = lmeDat, method = "ML", 
                     random =  list(rand1, rand2, rand4), 
                     correlation = corAR1())


#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", "W/O Year:Trophic Level", 
              "W/O Year:Initial S", "W/O Year:Study Duration",  "W/O Abundance", 
              "Year and abundance only", "Null model")

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_text, sort=TRUE, second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/richAIC reducedAbund.csv")

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: RICHNESS
### REDUCED DATASET, WITHOUT ABUNDANCE IN MODELS
############################################################
############################################################

lmeDat <- subDat
unique(lmeDat$Prediction)
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")
dim(lmeDat)

# Set up candidate model list
Cand.mod <- list()

# final full model 
Cand.mod[[1]] <- lme(fixed = richLN ~ 1 +  
	year0Z*Scale + year0Z*Prediction +
	year0Z*Trophic + year0Z*initialRichLN + 
	year0Z*durationZ, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())

# remove Scale interaction
Cand.mod[[2]] <- update(Cand.mod[[1]], ~. - year0Z:Scale)

# remove Prediction interaction 
Cand.mod[[3]] <- update(Cand.mod[[1]], ~. - year0Z:Prediction)

# remove Trophic interaction
Cand.mod[[4]] <- update(Cand.mod[[1]], ~. - year0Z:Trophic)

# remove initialRichLN interaction
Cand.mod[[5]] <- update(Cand.mod[[1]], ~. - year0Z:initialRichLN)

# remove duration interaction
Cand.mod[[6]] <- update(Cand.mod[[1]], ~. - year0Z:durationZ)

# year only
Cand.mod[[7]] <- lme(fixed = richLN ~ year0Z, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())

# null model
Cand.mod[[8]] <- lme(fixed = richLN ~ 1, 
	data = lmeDat, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
	
#create a vector of names to trace back models in set
mod_numbers <- paste("Cand.mod", 1:length(Cand.mod), sep=" ")	
mod_text <- c("Full model", "W/O Year:Scale", "W/O Year:Driver", "W/O Year:Trophic Level", 
              "W/O Year:Initial S", "W/O Year:Study Duration", "Year only", "Null model")

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_text, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
write.csv(mod.aicctab, "./output/richAIC reduced no abund.csv")

#generate AICc table
mod.aicctab <- aictab(cand.set= Cand.mod, modnames=mod_numbers, sort=TRUE, 
                      second.ord=FALSE) # second.ord =TRUE means AICc is used (not AIC)
print(mod.aicctab, digits=2, LL=TRUE)
