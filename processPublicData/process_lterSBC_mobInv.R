#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 150629
#################################################
# Dataset: Santa Barbara Coastal LTER
# www.sbc.lternet.edu
# Ongoing time-series of kelp forest community structure

# Quad and swath data
# Fixed species list
# Select mobile inverts only
# Aggregate small and large individuals (across quad and swath)

library(reshape2)
library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)

rm(list=ls(all=TRUE))

# Source functions
source("./divMetricF.R")

# Load SBC data
dat <- read.csv("./data/quad_swath_all_years_20140908.csv", 
                header=TRUE, na.strings= c("NA"))

dat$dateR <- as.Date(strftime(dat$DATE, format = "%Y-%m-%d"))
dat <- dat[with(dat, order(SITE, TRANSECT, dateR)), ]

# compare unique genus_species with codes
genSpF <- function(x) return(as.factor(paste(x$TAXON_GENUS, 
                                             x$TAXON_SPECIES, sep="_")))
dat$genSp <- genSpF(dat)

################################
# CREATE RAW DATA FILE 
################################

# Algae, sessile and mobile invertebrates
# 1m quads, and swaths
# Algal data not appropriate for biodiversity metrics, so remove

# -99999 refers to Macrocystis, therefore can just select inverts
dat2 <- droplevels(dat[dat$GROUP == "INVERT", ])
dim(dat2)
unique(dat2$SP_CODE)
unique(dat2$genSp)
# 55 species, 65 codes (because large and small individuals)

# Remove sessile species (better sampled as percent cover)
dat3 <- droplevels(filter(dat2, Mobility == "MOBILE"))

# REMOVE HAKA
dat4 <- droplevels(filter(dat3, SP_CODE != "HAKA"))
unique(dat4$SP_CODE)
unique(dat4$genSp)
summary(dat4$COUNT)
unique(dat4$COUNT)
glimpse(dat4)

# change -99999 to NA
COUNT2 <- dat4$COUNT
summary(COUNT2)
unique(COUNT2)
COUNT2[COUNT2 == -99999.0] <- NA
plot(dat4$COUNT ~ COUNT2, ylim = c(0, 350))

dat4$COUNT <- COUNT2

# now aggregate counts across quads and swaths
dat5 <- aggregate(COUNT ~ YEAR + MONTH + DATE + dateR +
                    SITE + TRANSECT + 
                    TAXON_GENUS + TAXON_SPECIES + genSp, 
                  data = dat4, sum)
summary(dat5)
dat5 <- dat5[with(dat5, order(SITE, TRANSECT, dateR)), ]
unique(dat5$COUNT)

unique(dat5$genSp)

# Need to fix Small Kelletia
dat6 <- filter(dat5, TAXON_GENUS != "Small Kelletia")

kellDat <- filter(dat5, TAXON_GENUS == "Small Kelletia")
head(kellDat)
kellDat$TAXON_GENUS <- "Kelletia"
kellDat$TAXON_SPECIES <- "kelletii"
kellDat$genSp <- "Kelletia_kelletii"

dat7 <- rbind(dat6, kellDat)
unique(dat7$genSp)

# now aggregate again
dat8 <- aggregate(COUNT ~ YEAR + 
                    SITE + TRANSECT + 
                    TAXON_GENUS + TAXON_SPECIES + genSp, 
                  data = dat7, sum)
unique(dat8$genSp)

# remove year 2000
dat9 <- filter(dat8, YEAR > 2000)
datL <- dat9

################################
# RAW DATASET IS READY, PREP FOR ANALYSIS
################################
# go from long to wide
names(datL)
datW <- dcast(datL, YEAR + SITE + TRANSECT ~ genSp, 
               value.var = "COUNT")
datW <- datW[with(datW, order(SITE, TRANSECT, YEAR)), ]

# remove rows with NAs - not sampled well
datW <- datW[complete.cases(datW), ]

# get rich, div, even, abund
library(vegan)
names(datW)
datW$rich <- specnumber(datW[, 4:36])
datW$div <- diversity(datW[, 4:36])
datW$even <- datW$div/(log(datW$rich))
datW$abund <- rowSums(datW[, 4:36])

# check for duplicates
styF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, x$YEAR, sep="_")))
datW$sty <- styF(datW)
which(duplicated(datW$sty))

# now look for replication within sites across year
stF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, sep="_")))
datW$st <- stF(datW)

# which site year combinations need to be removed so that I can calculate gamma div?
with(datW, table(st, YEAR))
# AQUE - remove AQUE_2, AQUE_6 (each missing one year)
# IVEE - remove IVEE 3-8
# SCTW - remove SCTW_1

unique(datW$st)
datW <- droplevels(filter(datW, st != "AQUE_2" &
                             st != "AQUE_6" &
                             st != "IVEE_3" &
                             st != "IVEE_4" &
                             st != "IVEE_5" &
                             st != "IVEE_6" &
                             st != "IVEE_7" &
                             st != "IVEE_8" &
                             st != "SCTW_1"))

# get site X year means for richness, diversity, evenness, and abundance; replicates
names(datW)
datW <- datW %>% rename(site = SITE, year = YEAR)

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
temp1 <- temp0[, 4:36] # MODIFY!
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
btsStudySub <- "SBCmobS"
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