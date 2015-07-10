#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Randomization procedure for diversity
# Author: Robin Elahi
# Date: 150629
#################################################

# General approach
# We now have 4 categories: negative, positive, neutral, or none.  The 'none' doesn't really count.  

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
head(divDat$year0Z)
tail(divDat$year0Z)

# diversity
dat <- divDat

# only one neutral prediction - lme won't run - need to use old classification
with(dat, table(Prediction))
dat$PredictionOld <- droplevels(as.factor(gsub("none", "neutral", dat$Prediction)))
unique(dat$PredictionOld)

# relevel predictionOld
dat$PredictionOld <- relevel(dat$PredictionOld, ref = "neutral")
with(dat, table(PredictionOld))

############################################################
############################################################
#Get p-value
getTrtP <- function(newDF){
  model.i <- lme(fixed = div ~ 1 +  
	year0Z*Scale + year0Z*PredictionOld +
	year0Z*Trophic + year0Z*initialDiv + 
	year0Z*durationZ, 
	data = newDF, method = "ML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1(), weights = varExp())

	
	overallP <- anova(model.i)[9,4] # overall p value
	levelP <- unname(summary(model.i)$tTable[11:12, 5]) 
	# p values for each level
	pVector <- c(overallP, levelP)
	return(pVector)
}

############################################################
############################################################
# Function to switch negative and positive to neutral, neutral to either
# create vector of new predictions corresponding to length of old prediction
switchFactor3 <- function(x) {
	predictionOriginal <- unique(x)
	if(predictionOriginal == "negative") {predictionNew <- "neutral"}
	if(predictionOriginal == "positive") {predictionNew <- "neutral"}
	if(predictionOriginal == "neutral") 
		{predictionNew <- sample(c("negative", "positive"), 1)}
	predictionVector <- rep(predictionNew, length(x))
	return(predictionVector)}

############################################################
############################################################
# Function to change more than 1 time series' predictions
# wrap into function
tsNames <- unique(dat[dat$Prediction != "none", ]$subSiteID)
tsNames <- unique(dat$subSiteID)
tsNames
tsN <- length(tsNames); tsN

switcharoo <- function(dat, newPredictions) {
	randSample <- sample(tsNames, newPredictions)

	# create subset to change
	randChange <- dat[dat$subSiteID %in% randSample, ]
	
	# create subset to keep
	randKeep <- dat[!dat$subSiteID %in% randSample, ]

	# get list of new predictions
	predictionList <- tapply(X = randChange$PredictionOld, 
	INDEX = randChange$subSiteID, FUN = switchFactor3)

	#need to change from list to a vector/data frame
	predictionNew <- factor(matrix(unlist(predictionList)))

	# create new changed subset of data
	randChange2 <- randChange
	randChange2$PredictionOld <- predictionNew

	# create the new data frame
	newDF <- rbind(randKeep,  randChange2)

	# run the model and get the pValue
	getTrtP(newDF)
}

# switcharoo(dat, 4)

############################################################
############################################################
# now run this n times for 1:13
predictionChanges <- c(1:13)
predictionChanges

AllReps <- predictionChanges
N <- length(AllReps); N

pValsList <- vector("list", length = N)
pValsList
names(pValsList) <- factor(predictionChanges)
str(pValsList)

pValsN <- 10 # number of times to replicate

# estimated time
unitTime <- 14 # 67 seconds per lme call
unitTime/60 * pValsN * N # in minutes 

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
names(pValsDF) <- c("predictionChanges", "overall", "negative", "positive")
pValsDF$nChanges <- as.numeric(pValsDF$predictionChanges)
pValsDF
pValsLong <- melt(pValsDF, id.vars = c("predictionChanges", "nChanges"))
pValsLong

write.csv(pValsLong, "./output/pValsLongDiversity_150620.csv")
