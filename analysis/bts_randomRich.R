#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Randomization procedure for richness
# Author: Robin Elahi
# Date: 150629
#################################################

# General approach
# We now have 4 categories: negative, positive, neutral, or none. 
# The 'none' doesn't really count.

# For each iteration, change a single time series within each category as follows:
# Negative --> neutral
# Positive --> neutral
# Neutral --> positive or negative (coin flip)
# Then re-assess model fit, and note the proportion of times (out of 100?) the effect of driver is no longer 'significant'

# I could repeat the following, but instead change two at once, three at once, all the way up to the total # of time series that were categorized.  

# Then I could plot: proportion of times driver is no longer significant, against the # of 'switches'.  Then we'd have an estimate of the tolerance to categorizing error.  


rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# get data for analysis
source("./R/bts_dataPrep.R")

# detach dplyr
detach("package:dplyr", unload = TRUE)

library(nlme)
library(reshape2)
library(plyr)

# Random effects for LME models
source("./R/bts_lme_models.R")

# make sure the data are ordered by time within time-series (
# necessary for autocorrelation 
head(richDat$year0Z)
tail(richDat$year0Z)

# richness
dat <- richDat
dat$Prediction <- relevel(dat$Prediction, ref = "none")






############################################################
############################################################
#Write function to get p-value from LME call

getTrtP <- function(newDF){
  model.i <- lme(fixed = richLN ~ 1 +  
	year0Z*Scale + year0Z*Prediction +
	year0Z*Trophic + year0Z*initialRichLN + 
	year0Z*durationZ, 
	data = newDF, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
	
	overallP <- anova(model.i)[9,4] # overall p value
	levelP <- unname(summary(model.i)$tTable[12:14, 5]) # p values for each level
	pVector <- c(overallP, levelP)
	return(pVector)
}

############################################################
############################################################
# Function to switch negative and positive to neutral, neutral to either
# create vector of new predictions corresponding to length of old prediction
switchFactor3 <- function(x) {
	predictionOld <- unique(x)
	if(predictionOld == "negative") {predictionNew <- "neutral"}
	if(predictionOld == "positive") {predictionNew <- "neutral"}
	if(predictionOld == "neutral") 
		{predictionNew <- sample(c("negative", "positive"), 1)}
	predictionVector <- rep(predictionNew, length(x))
	return(predictionVector)}
		
############################################################
############################################################
# Function to change more than 1 time series' predictions
# wrap into function
tsNames <- unique(dat[dat$Prediction != "none", ]$subSiteID)
tsNames
tsN <- length(tsNames); tsN

switcharoo <- function(dat, newPredictions) {
	randSample <- sample(tsNames, newPredictions)

	# create subset to change
	randChange <- dat[dat$subSiteID %in% randSample, ]
	
	# create subset to keep
	randKeep <- dat[!dat$subSiteID %in% randSample, ]

	# get list of new predictions
	predictionList <- tapply(X = randChange$Prediction, 
	INDEX = randChange$subSiteID, FUN = switchFactor3)

	#need to change from list to a vector/data frame
	predictionNew <- factor(matrix(unlist(predictionList)))

	# create new changed subset of data
	randChange2 <- randChange
	randChange2$Prediction <- predictionNew

	# create the new data frame
	newDF <- rbind(randKeep,  randChange2)

	# run the model and get the pValue
	getTrtP(newDF)
}

############################################################
############################################################
# run for 36 alternative classifications 
predictionChanges <- seq(1, 72, by = 2)

AllReps <- predictionChanges
N <- length(AllReps); N

pValsList <- vector("list", length = N)
pValsList
names(pValsList) <- factor(predictionChanges)
str(pValsList)

pValsN <- 10 # number of times to replicate

# estimated time
unitTime <- 10 # 7 seconds per lme call
unitTime/60 * pValsN * N 

AllReps[2]

# use trycatch
for (i in 1:N) {
  tryCatch({
  predictionChange <- AllReps[i]
  pVals <- replicate(n = pValsN, switcharoo(dat, 
                                            newPredictions = predictionChange))
  pValsList[[i]] <- unname(t(pVals))
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
}


pValsList

# RESHAPE DATA
pValsDF <- ldply(pValsList)
names(pValsDF) <- c("predictionChanges", "overall", "negative", "neutral", "positive")
pValsDF$nChanges <- as.numeric(pValsDF$predictionChanges)
pValsDF
pValsLong <- melt(pValsDF, id.vars = c("predictionChanges", "nChanges"))
pValsLong

write.csv(pValsLong, "./output/pValsLong10reps_150620.csv")

########################
# MISSING '65' values
pVals65 <- replicate(n = pValsN, switcharoo(dat, 
                                            newPredictions = 65))
pVals65
pValsList65 <- unname(t(pVals65))
pValsList65

# RESHAPE DATA
pValsDF <- as.data.frame(pValsList65)
pValsDF
names(pValsDF) <- c("overall", "negative", "neutral", "positive")

nChanges <- rep(65, 10)
predictionChanges <- rep(65, 10)

pValsDF2 <- cbind(predictionChanges, nChanges, pValsDF)
pValsDF2

pValsLong <- melt(pValsDF2, id.vars = c("predictionChanges", "nChanges"))
pValsLong

write.csv(pValsLong, "./output/pValsLong10reps_65_150620.csv")
