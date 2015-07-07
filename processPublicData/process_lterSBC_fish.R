#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 150629
#################################################
# Dataset: Santa Barbara Coastal LTER
# www.sbc.lternet.edu
# Ongoing time-series of kelp forest community structure

# Abundance and size of reef fish
# Aggregate cryptic fish and fish on each transect
# Only use autumn data
# Lump species they had lumped in 2005

library(reshape2)
library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)

rm(list=ls(all=TRUE))

# Source functions
source("./divMetricF.R")

# Load SBC data from lter website
dat <- read.csv("./data/all_fish_all_years_20140903.csv", 
                header=TRUE, na.strings= c("NA", "NULL"))
dim(dat)
# REMOVE OTRI - SP_CODE
dat <- filter(dat, SP_CODE != "OTRI")
# Select autumn data only
dat <- filter(dat, SURVEY_TIMING == "A")

# Dates	
dat$dateR <- as.Date(strftime(dat$date, format = "%Y-%m-%d"))
dat <- dat[with(dat, order(SP_CODE, site, transect, dateR)), ]

# remove data with count =  -99999, 5 observations
dat <- filter(dat, COUNT != -99999)

################################
# CREATE RAW DATA FILE 
################################
# The raw data have sizes, but I want counts
dat2 <- aggregate(COUNT ~ year + month + date + dateR + 
                    SURVEY_TIMING  + 
                    site + transect + 
                    + TAXON_GENUS + TAXON_SPECIES + SP_CODE, 
                  data = dat, sum)

dat2 <- droplevels(dat2[with(dat2, order(SP_CODE, site, transect, dateR)), ])
summary(dat2)
################################
# Remove row IF genus is no fish,
dat3 <- filter(dat2, TAXON_GENUS != "No fish")
dim(dat3)

# create temp data frame to extract sites with no fish
tempDat <- filter(dat2, TAXON_GENUS == "No fish")
dim(tempDat) # only one transect had no fish 	
tempDat <- filter(tempDat, COUNT == 1)

################################
# If TAXON_GENUS = -99999, then remove row
# Must do this before adding the one observation of 'no fish'
dat4 <- dat3[dat3$TAXON_GENUS != "-99999", ]
dim(dat4)

# Add no fish back in
dat5 <- rbind(dat4, tempDat)			

styF <- function(x) return(as.factor(paste(x$site, x$transect, x$year, sep="_")))
dat5$sty <- styF(dat5)

# Check mohawk1 in 2009
mohk1 <- filter(dat5, site == "MOHK" & transect == "1" & year == "2009")
mohk1 <- filter(dat5, sty == "MOHK_1_2009")
mohk2 <- filter(mohk1, date == "2009-07-23")
mohk2 <- slice(mohk1, c(1:8, 12))
mohk2 # good data for mohawk1 in 2009
# remove original mohawk 1 2009, then add back in the revised version
dat6 <- filter(dat5, sty != "MOHK_1_2009")
dim(dat5); dim(dat6)
dat7 <- rbind(dat6, mohk2)
dim(dat7)

################################
### Abundance
dat7$Abundance <- dat7$COUNT
dat7 <- dat7[with(dat7, order(site, transect, dateR)), ]

siteTranF <- function(x) return(as.factor(paste(x$site, x$transect, sep="_")))
dat7$siteTran <- siteTranF(dat7)
unique(dat7$siteTran) 

################################
# lump species that were lumped < 2005
genSpF <- function(x) return(as.factor(paste(x$TAXON_GENUS, x$TAXON_SPECIES, sep="_")))

dat7$genSp <- genSpF(dat7)

dat04 <- droplevels(filter(dat7, year < 2005))
dat05 <- droplevels(filter(dat7, year > 2004))
dim(dat05)

spp04 <- unique(dat04$genSp)
spp05 <- unique(dat05$genSp)

spp04 <- paste(c(spp04[order(spp04)], rep('null', 20)))
spp05 <- paste(spp05[order(spp05)])

tempDat <- data.frame(spp04, spp05)

# lump Sebastes flavidus, miniatus, paucispinis, rastrelliger, 
# and serriceps to spp
dat8 <- dat7
dat8$genSp <- gsub('Sebastes_flavidus', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_miniatus', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_paucispinis','Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_rastrelliger', 'Sebastes_spp.', dat8$genSp)
dat8$genSp <- gsub('Sebastes_serriceps', 'Sebastes_spp.', dat8$genSp)
unique(dat8$genSp)
unique(dat7$genSp)

# now create new species column from new genSp column
# This is a bit of a hack
species2 <- unlist(strsplit(dat8$genSp, "_"))
length(species2)
# here i take the even values only
species3 <- species2[seq(2, length(species2), by = 2)] 

dat9 <- dat8
dat9$TAXON_SPECIES <- species3
unique(dat9$genSp)

################################
# RAW DATASET IS READY, PREP FOR ANALYSIS
################################
# go from long to wide
datL <- dat9

# rename columns to match
datL <- rename(datL, YEAR = year, SITE = site, TRANSECT = transect)
names(datL)

datW <- dcast(datL, YEAR + SITE + TRANSECT ~ genSp, 
              value.var = "Abundance")

names(datW)
summary(datW)
datW <- datW[with(datW, order(SITE, TRANSECT, YEAR)), ]

# remove NO FISH column
datW <- datW[, -33]
names(datW)

# get rich, div, even, abund
library(vegan)
names(datW)
datW$rich <- specnumber(datW[, 4:58])
datW$div <- diversity(datW[, 4:58])
datW$even <- datW$div/(log(datW$rich))
datW$abund <- rowSums(datW[, 4:58])

# check for duplicates
styF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, x$YEAR, sep="_")))
datW$sty <- styF(datW)
which(duplicated(datW$sty))

# now look for replication within sites across year
stF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, sep="_")))
datW$st <- stF(datW)

# remove year 2000
datW <- filter(datW, YEAR > 2000)

# which site year combinations need to be removed so that I can calculate gamma div?
with(datW, table(st, YEAR))

unique(datW$st)
datW2 <- droplevels(filter(datW, 
                            st != "AQUE_6" &
                            st != "BULL_2" &
                            st != "BULL_4" &
                            st != "BULL_5" &
                            st != "BULL_7" &
                            st != "BULL_8" &
                            st != "CARP_2" &
                            st != "CARP_4" &
                            st != "CARP_8" &
                            st != "IVEE_3" &
                            st != "IVEE_4" &
                            st != "IVEE_5" &
                            st != "IVEE_6" &
                            st != "IVEE_7" &
                            st != "IVEE_8" &
                            st != "SCDI_1"))

with(datW2, table(st, YEAR))

# select good ahnd, golb, bull
ahnd <- droplevels(filter(datW2, SITE == "AHND" &
                             YEAR > 2006))
golb <- droplevels(filter(datW2, SITE == "GOLB" &
                            YEAR > 2005))
bull <- droplevels(filter(datW2, SITE == "BULL" &
                            YEAR < 2011))

# remove ahdn, golb, bull
datW3 <- filter(datW2, SITE != "AHND" & 
                  SITE != "GOLB" &
                  SITE != "BULL")

datW4 <- rbind(datW3, ahnd, golb, bull)
dim(with(datW4, table(st, YEAR)))

# RENAME
datW <- datW4 %>% rename(site = SITE, year = YEAR)

# get site X year means for richness, diversity, evenness, and abundance; replicates
names(datW)
sdF <- function(x) return(as.factor(paste(x$site, x$year, sep="_")))
datW$sd <- sdF(datW)

# dplyr
alphDF <- datW %>% group_by(site, year, sd) %>% 
  summarise(rich = mean(rich), div = mean(div), even = mean(even), 
            abund = mean(abund), nReps = n()) %>% 
  ungroup()
alphDF

### Calculate sums per site, to get gamma scale metrics
names(datW)

# Create a uniform data frame
temp0 <- datW # MODIFY
temp1 <- temp0[, 4:58] # MODIFY!
temp2 <- select(temp0, site, year, sd)
temp3 <- cbind(temp1, temp2)
head(temp3)

gamDF <- temp3 %>% group_by(site, year, sd) %>% 
  summarise_each(funs(sum)) %>% ungroup()
head(gamDF)

# Gamma - rich, div, even, abund
gamDF <- divMetF(gamDF, 4, dim(gamDF)[2])
gamDF$rich

# create dataframes

# study traits
btsStudyName <- "SBC"
btsStudySub <- "SBCfish"
alphaSQM <- 80
gammaSQM <- NA
abundType <- "Abundance"

# ALPHA
btsDFalpha <- data.frame(
  studyname = btsStudyName, 
  studySub = btsStudySub, 
  scale = "plot", 
  sppCode = "R", 
  siteScale = "alpha", 
  site = alphDF$site, 
  subSiteID = "", 
  nest1 = NA, 
  nest2 = NA, 
  date = as.Date(paste(alphDF$year, 08, 01, sep = "-")), 
  date.no = as.Date(paste(alphDF$year, 08, 01, sep = "-")),
  rich = alphDF$rich, 
  div = alphDF$div, 
  even = alphDF$even, 
  rep = alphDF$site, 
  remove = 0, 
  removeNote = NA, 
  abund = alphDF$abund, 
  abundSQM = alphDF$abund/alphaSQM,
  sqmOriginal = alphaSQM, 
  nReps = alphDF$nReps[1]
)

head(btsDFalpha)
btsDFalpha$subSiteID <- with(btsDFalpha, paste(studySub, site, sppCode, sep = "_"))

# GAMMA
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
  date = as.Date(paste(alphDF$year, 08, 01, sep = "-")), 
  date.no = as.Date(paste(alphDF$year, 08, 01, sep = "-")),
  rich = gamDF$rich, 
  div = gamDF$div, 
  even = gamDF$even, 
  rep = gamDF$site, 
  remove = 0, 
  removeNote = NA, 
  abund = gamDF$abund, 
  abundSQM = gamDF$abund/(alphaSQM*alphDF$nReps[1]),
  sqmOriginal = alphaSQM*alphDF$nReps[1], 
  nReps = alphDF$nReps[1]
)

head(btsDFgamma)
btsDFgamma$subSiteID <- with(btsDFgamma, paste(studySub, site, sppCode, sep = "_"))

btsDF <- rbind(btsDFalpha, btsDFgamma)
head(btsDF)
tail(btsDF)







