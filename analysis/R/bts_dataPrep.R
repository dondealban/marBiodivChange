#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Creating data subsets for the analyses found in:
# Author: Robin Elahi
# Date: 150629
#################################################

library(dplyr)

#Load in the master data set
fullDat <- read.csv("./data/TableS4.csv", header=TRUE, na.strings="NA")

with(fullDat, c(length(unique(studyName)), length(unique(studySub)), 
                length(unique(subSiteID)), length(studyName)))
unique(fullDat$site)
glimpse(fullDat)

# change trophic level to a factor (from integer)
fullDat$Trophic <- as.factor(fullDat$Trophic)
with(fullDat, c(length(unique(studyName)), length(unique(studySub)), 
                length(unique(subSiteID)), length(studyName), length(unique(site))))

# order the data by time
fullDat <- fullDat[with(fullDat, order(studyName, studySub, subSiteID, dateR)), ]
head(fullDat)

### Create subsets of the full dataset for analysis

# complete richness dataset
richDat <- droplevels(fullDat[complete.cases(fullDat$rich), ])
dim(richDat)
names(richDat)
length(unique(richDat[, 2]))
# get # of studies, substudies, and time series
with(richDat, c(length(unique(studyName)), length(unique(studySub)), 
                length(unique(subSiteID)), length(studyName)))

# now create richness dataset that also has diversity
richDivDat <- droplevels(richDat[complete.cases(richDat$div), ])
dim(richDivDat)
with(richDivDat, c(length(unique(studyName)), length(unique(studySub)), 
               length(unique(subSiteID)), length(studyName)))

# now create richness dataset that also has abundance
subDat <- droplevels(richDat[complete.cases(richDat$abundFinal), ])
dim(subDat)

# get # of studies, substudies, and time series
with(subDat, c(length(unique(studyName)), length(unique(studySub)), 
               length(unique(subSiteID)), length(studyName)))

####################################################
# create diversity dataset
divDat <- droplevels(fullDat[complete.cases(fullDat$div), ])
dim(divDat)
# get # of studies, substudies, and time series
with(divDat, c(length(unique(studyName)), length(unique(studySub)), 
               length(unique(subSiteID)), length(studyName)))


# now create diversity dataset that also has abundance
divDatSub <- droplevels(divDat[complete.cases(divDat$abundFinal), ])
dim(divDatSub)

# get # of studies, substudies, and time series
with(divDatSub, c(length(unique(studyName)), length(unique(studySub)), 
                  length(unique(subSiteID)), length(studyName)))



