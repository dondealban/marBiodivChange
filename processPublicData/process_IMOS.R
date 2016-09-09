#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 150629
#################################################
# Dataset: IMOS
# https://imos.aodn.org.au/imos123/
# IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance

# PH4 was sampled with smaller net size (20cm, as opposed to 60cm), 
# and sampled on the way up and down (as opposed to only down)
# Adult copepods only
# Aggregate abundance at genus level because there were many
# taxonomic adjustments at the species level throughout the time-series

##' Data were downloaded from IMOS website on 28 May 2015

library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)

rm(list=ls(all=TRUE))

# Source functions
source("./divMetricF.R")

#Load in the master data set
dat <- read.csv("./imosDownload_150528/IMOS_National_Reference_Station_(NRS)_-_Zooplankton_Abundance.csv", 
                header = TRUE, na.strings = "", skip = 54)
head(dat)
str(dat) # 10 stations
names(dat)
summary(dat)
dim(dat)
with(dat, table(TAXON_GROUP))

# create row ID
dat$ID <- seq(1:dim(dat)[1])

unique(dat$GENUS)
genus2 <- as.factor(sub(" +$", "", dat$GENUS)) # trim trailing white space
unique(genus2)

unique(dat$SPECIES)
species2 <- as.factor(sub(" +$", "", dat$SPECIES)) # trim trailing white space
unique(species2)

genSpF<-function(x) return(as.factor(paste(x$GENUS, x$SPECIES, sep="_")))
dat$genSp <- genSpF(dat)
unique(dat$genSp) 
unique(dat$genus_and_species)

# Focus on adults due to identification problems
adultDat <- droplevels(dat[dat$LIFE_STAGE == "ADULT", ])
dim(adultDat)
unique(adultDat$LIFE_STAGE) 
unique(adultDat$TAXON_NAME)
unique(adultDat$taxon_group) 
unique(dat$genSp) 

with(adultDat, table(TAXON_GROUP))
names(adultDat)
adultDat <- adultDat[with(adultDat, order(genSp, STATION_NAME, UTC_TRIP_START_TIME)), ]

unique(adultDat$GENUS)

# AGGREGATE TO GENUS LEVEL
names(adultDat)
adultDat2 <- aggregate(TAXON_PER_M3 ~ STATION_NAME +
                        UTC_TRIP_START_TIME + 
                        LONGITUDE + LATITUDE + GENUS, 
                      data = adultDat, sum)

################################
# RAW DATASET IS READY, PREP FOR ANALYSIS
################################
# go from long to wide
datL <- adultDat2
names(datL)

# rename columns to match
datL <- rename(datL, DATE = UTC_TRIP_START_TIME, 
               SITE = STATION_NAME, COUNT = TAXON_PER_M3)
names(datL)
head(datL)

library(reshape2)
datW <- dcast(datL, DATE + SITE ~ GENUS, 
              value.var = "COUNT")
names(datW)
summary(datW)
datW <- datW[with(datW, order(SITE, DATE)), ]

# now I need to replace NA with 0 in abundance columns
head(datL) # NAs are incorporated during dcast
datW[is.na(datW)] <- 0
head(datW)

# Gamma - rich, div, even, abund
names(datW)
datW <- divMetF(datW, 3, 101)

# DEAL WITH TIME
date1 <- substr(datW$DATE, 1, 10)
str(date1)
datW$dateR <- as.Date(strptime(date1, format = "%Y-%m-%d"))

# check for duplicates
sdF <- function(x) return(as.factor(paste(x$SITE, x$dateR, sep="_")))
datW$sd <- sdF(datW)
which(duplicated(datW$sty))

# REMOVE PH4 and SOTS
unique(datW$SITE)

datW <- filter(datW, SITE != "Southern Ocean Time Series" &
                 SITE != "Port Hacking 4")

qplot(dateR, rich, data = datW, 
      facets = ~SITE, geom = "line")

gamDF <- datW

# study traits
btsStudyName <- "IMOS"
btsStudySub <- "IMOS"
alphaSQM <- NA
gammaSQM <- 0.282743339

# GAMMA
gamDF <- gamDF %>% rename(site = SITE)

btsDFgamma <- data.frame(
  studyname = btsStudyName, 
  studySub = btsStudySub, 
  scale = "siteObs", 
  sppCode = "S", 
  siteScale = "gamma", 
  site = gamDF$site, 
  subSiteID = "", 
  nest1 = NA, 
  nest2 = NA, 
  date = gamDF$dateR, 
  date.no = gamDF$dateR, 
  rich = gamDF$rich, 
  div = gamDF$div, 
  even = gamDF$even, 
  rep = gamDF$site, 
  remove = 0, 
  removeNote = NA, 
  abund = gamDF$abund, 
  abundSQM = NA, 
  sqmOriginal = gammaSQM, 
  nReps = NA
)

head(btsDFgamma)
btsDFgamma$subSiteID <- with(btsDFgamma, paste(studySub, site, sppCode, sep = "_"))

btsDF <- rbind(btsDFgamma)
head(btsDF)
tail(btsDF)
