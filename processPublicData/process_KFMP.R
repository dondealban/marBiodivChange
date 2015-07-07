#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 150629
#################################################
# Dataset: KFMP
# http://www.esajournals.org/doi/abs/10.1890/13-0562R.1
# KFMP (kelp forest monitoring program - channel islands)
# RDFC (roving diver fish counts)
# all fish species counted on a 30 minute dive along a site transect
# multiple 'experts' sometimes surveyed the same site transect for a given time point

library(reshape2)
library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)

rm(list=ls(all=TRUE))

# Source functions
source("./divMetricF.R")

rm(list=ls(all=TRUE)) 

#Load in the master data set
dat <- read.csv("./data/RDFC data_150529.csv", header = TRUE,  na.strings = "NA")
dat$dateR<-as.Date(dat$date.no, origin="1904-01-01")
head(dat)

##########################################################

unique(dat$ScientificName) # 164
# need to collapse 'all' with 'juvenile', 'adult', 'female', 'male'

# use gsub to remove f and m
# use trim to get rid of trailing spaces
# tapply to sum the abundance of each species
# then dcast

taxon <- gsub(" all|juvenile|adult|male|female", "", dat$ScientificName)
unique(taxon)
taxon <- sub(" +$", "", taxon) # trim trailing white space
taxon <- sub(",+$", "", taxon) # trim trailing commas
taxon <- sub(" +$", "", taxon) # trim trailing white space

taxon <- as.factor(taxon)
taxaNames <- unique(taxon) # 109 unique
taxaNames[order(taxaNames)]

##########################################################

dat1 <- cbind(dat, taxon)
names(dat1)

# drop baitfish and larval fish
dat2 <- dat1[dat1$taxon != "baitfish unidentified" &
						dat1$taxon != "larval fish spp.", ]
dim(dat2); dim(dat1)

# also need to lump Bathymasteridae and spp.
practice <- unique(dat2$taxon)
practice
grepl("Bathymasteridae", practice)
practice2 <- ifelse(grepl("Bathymasteridae", practice) == TRUE, 
									"crap", as.character(practice))
practice2

# now do this for the dataset
taxon2 <- ifelse(grepl("Bathymasteridae", dat2$taxon) == TRUE, 
									"Bathymasteridae", as.character(dat2$taxon))
length(taxon2)
dat2$taxon2 <- as.factor(taxon2)
str(dat2$taxon2)

# Counts are the number observed per 30 minutes - NOT DENSITY
# And, there were some cases in more recent years 
# when counts were not possible
# THEREFORE: PRESENCE/ABSENCE DATA
datLong <- dat2
head(datLong)
unique(datLong$Abundance)
summary(datLong)

abund1 <- as.character(datLong$Abundance)
summary(abund1)
abund1
str(abund1)
abund1[is.na(abund1)] <- 0

datLong$PA <- ifelse(abund1 == "0", 0, 1)
tail(datLong)

###
# lump Sebastes with '/'
unique(datLong$taxon2)

datL <- datLong

datL$taxon2 <- gsub('Sebastes serranoides/flavidus', 'Sebastes spp.', 
                   datL$taxon2)
datL$taxon2 <- gsub('Sebastes atrovirens/carnatus/caurinus/chrysomelas', 
                    'Sebastes spp.', 
                    datL$taxon2)
datL$taxon2 <- gsub('Sebastes chrysomelas/carnatus', 
                    'Sebastes spp.', 
                    datL$taxon2)

levels(as.factor(datL$taxon2))
unique(datLong$taxon2)
##########################################################
# need to remove duplicate taxa for each site x time x expert group
datL$dummy <- rep(1, dim(datL)[1])
head(datL)

datL2 <- aggregate(PA ~ Site + Year + dateR + ObserverNumber + 
                        taxon2 + dummy, data = datL, sum)

unique(datL2$taxon2) # correct # of taxa
sum(datL2$dummy) # the correct # of observations
head(datL2)

################################
# RAW DATASET IS READY, PREP FOR ANALYSIS
################################
# now i can melt the data
datWide <- dcast(datL2, Site + Year + dateR + ObserverNumber ~ taxon2, value.var = "PA")
summary(datWide)
dim(datWide)
names(datWide)

##########################################################
# check for duplicates
sdoF <- function(x) return(as.factor(paste(x$Site, x$dateR, x$ObserverNumber, 
                                           sep="_")))
datWide$sdo <- sdoF(datWide)
head(datWide)
which(duplicated(datWide$sdo))
summary(datWide)

# now I need to replace values greater than 0 with 1 in abundance columns
# and NAs with 0
datWide[is.na(datWide)] <- 0
summary(datWide)
names(datWide)

sppMat <- datWide[, 5:107]
sppMat[1:100, 19]
sppMat2 <- ifelse(sppMat == 0, 0, 1)
sppMat2[1:100, 19]
summary(sppMat2)

datFinal <- cbind(datWide[, 1:4], sppMat2)
head(datFinal)
dim(datFinal)
datFinal$rich <- specnumber(datFinal[, 5:107])

# take the average richness per site and date combo
datFinal2 <- datFinal %>% group_by(Site, dateR) %>%
  summarise(meanR = mean(rich))

qplot(dateR, rich, data = datFinal, 
      facets = ~Site, geom = "point")
qplot(dateR, meanR, data = datFinal2, 
      facets = ~Site, geom = "point")

with(datFinal, table(Year, ObserverNumber))
names(datFinal2)
summary(datFinal2)

##############################
gamDF <- datFinal2
# study traits
btsStudyName <- "KFMP"
btsStudySub <- "KFMP"
alphaSQM <- NA
gammaSQM <- 2000

# GAMMA
gamDF <- gamDF %>% rename(site = Site, rich = meanR)

btsDFgamma <- data.frame(
  studyname = btsStudyName, 
  studySub = btsStudySub, 
  scale = "siteObs", 
  sppCode = "S", 
  siteScale = "gamma", 
  site = paste("KFMP", gamDF$site, sep = ""), 
  subSiteID = "", 
  nest1 = NA, 
  nest2 = NA, 
  date = gamDF$dateR, 
  date.no = gamDF$dateR, 
  rich = gamDF$rich, 
  div = NA, 
  even = NA, 
  rep = gamDF$site, 
  remove = 0, 
  removeNote = NA, 
  abund = NA, 
  abundSQM = NA, 
  sqmOriginal = gammaSQM, 
  nReps = NA
)

head(btsDFgamma)
btsDFgamma$subSiteID <- with(btsDFgamma, paste(studySub, site, sppCode, sep = "_"))

btsDF <- rbind(btsDFgamma)
head(btsDF)
tail(btsDF)

