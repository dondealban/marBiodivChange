#################################################
# Prepare raw data for biodiversity change database
# Author: Robin Elahi
# Date: 150629
#################################################
# Dataset: Santa Barbara Coastal LTER
# www.sbc.lternet.edu
# Ongoing time-series of kelp forest community structure
# Using updated UPC data from Shannon Harrer
# Shannon Harrer <harrer@ucsb.edu>

# Point contact of sessile algae and invertebrates

library(reshape2)
library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)
library(vegan)

rm(list=ls(all=TRUE))

# Source functions
source("./divMetricF.R")

# Load SBC data
dat <- read.csv("./data/cover_all_years_20150618.csv", header=TRUE) 
# can't use NA b/c species code

dat$dateR <- as.Date(strftime(dat$DATE, format = "%Y-%m-%d"))
dat <- dat[with(dat, order(SITE, TRANSECT, dateR)), ]

# Create genSp
genSpF<-function(x) return(as.factor(paste(x$TAXON_GENUS, x$TAXON_SPECIES, sep="_")))
dat$genSp <- genSpF(dat)
unique(dat$genSp)

# REMOVE DMH - dead macro holdfast
dat2 <- filter(dat, SP_CODE != "DMH")
# REMOVE BLD - kelp juveniles
dat2 <- filter(dat2, SP_CODE != "BLD")

# REMOVE SUBSTRATE CODES
names(dat2)
unique(dat2$GROUP)
dat3 <- dat2[dat2$GROUP != "SUBSTRATE", ]
head(dat3)
unique(dat3$genSp)

dat3 <- dat3[with(dat3, order(TAXON_GENUS, TAXON_SPECIES, SITE, TRANSECT, dateR)), ]


################################
# NEED TO ACCOUNT FOR TAXONOMIC
# CHANGES OVER TIME
################################


genSp <- dat3 %>% group_by(SP_CODE, TAXON_GENUS, TAXON_SPECIES) %>%
  summarise(sum = sum(PERCENT_COVER))
genSp
length(unique(genSp$SP_CODE))
write.csv(genSp, 'genSp.csv')

unique(dat3$NOTES)
noteDat <-  droplevels(dat3[dat3$NOTES != "-99999", ])
summary(noteDat)
unique(noteDat$SP_CODE)

# concatenate species and note
speciesNote <- with(noteDat, paste(SP_CODE, NOTES, TAXON_GENUS, 
                                   TAXON_SPECIES, sep = "_"))
unique(speciesNote)

# SO - FIRST STEP IS TO RENAME THESE APPROPRIATELY
# AND THEN SUM

code1 <- dat3$SP_CODE

code2 <- gsub("CRYP", replace = "BF", code1)
code3 <- gsub("EUPO", replace = "SABW", code2)
code4 <- gsub("RP", replace = "R", code3)
code5 <- gsub("STIN", replace = "R", code4)
code6 <- gsub("ANSO", "ANSP", code5)
code7 <- gsub("DIAT", "FB", code6)
code8 <- gsub("EUCL", "UT", code7)
code9 <- gsub("GEL", "GR", code8)
code10 <- gsub("GYSP", "R", code9)
c11 <- gsub("OBSP", "UIH", code10)
c12 <- gsub("PHSP", "UNAN", c11)
c13 <- gsub("PLUM", "UIH", c12)
c14 <- gsub("PRSP", "R", c13)
c15 <- gsub("PYST", "UT", c14)

unique(c15)
length(unique(c15))

dat3$SP2 <- c15

# REMOVE FTHR, CYOS
dat4 <- filter(dat3, SP2 != "FTHR" & SP2!= "CYOS")

# NEXT STEP, ADD PERCENT COVER FOR SAME SP-CODES
unique(dat4$GROUP)
summary(dat4$PERCENT_COVER)
unique(dat4$PERCENT_COVER)
dat4$pc <- as.numeric(gsub("-99999", "0", dat4$PERCENT_COVER))
unique(dat4$pc)
summary(dat$pc)
summary(dat4)

dat5 <- aggregate(pc ~ YEAR + MONTH + dateR + 
                    SITE + TRANSECT + GROUP + 
                    SP2, data = dat4, sum)

length(unique(dat5$SP2))

dat6 <- dat5[with(dat5, order(SITE, TRANSECT, dateR, SP2)), ]

# remove 2000
datL <- filter(dat6, YEAR > 2000)

# separate inverts and algae
invert <- datL %>% filter(GROUP == "INVERT")
algae <- datL %>% filter(GROUP != "INVERT")

################################
# RAW DATASET IS READY, PREP FOR ANALYSIS
################################

glimpse(invert)
glimpse(algae)

# go from long to wide
invW <- dcast(invert, SITE + TRANSECT + dateR + YEAR ~ SP2, 
              value.var = "pc", sum)
head(invW)

algW <- dcast(algae, SITE + TRANSECT + dateR + YEAR ~ SP2, 
              value.var = "pc", sum)
head(algW)

### Algae ###
# now look for replication within sites across year
stF <- function(x) return(as.factor(paste(x$SITE, x$TRANSECT, sep="_")))
algW$st <- stF(algW)

# which site year combinations need to be removed so that I can calculate gamma div?
with(algW, table(st, YEAR))

# some year x st combos have two - these are mistakes
dubs <- algW %>% group_by(YEAR, st) %>%
  summarise(n = n()) %>% ungroup() %>%
  filter(n == 2)
dubs

styF <- function(x) return(as.factor(paste(x$st, x$YEAR, sep="_")))
dubs$sty <- styF(dubs)
dubs

algW$sty <- styF(algW)

# Now I have a matching column
tempDat <- filter(algW, sty == "MOHK_2_2002" |
                    sty == "ABUR_1_2003" |
                    sty == "BULL_1_2003" |
                    sty == "GOLB_2_2003" |
                    sty == "IVEE_1_2003")

algW2 <- filter(algW, sty != "MOHK_2_2002" &
                  sty != "ABUR_1_2003" &
                  sty != "BULL_1_2003" &
                  sty != "GOLB_2_2003" &
                  sty != "IVEE_1_2003")
tempDat
tempDat2 <- tempDat[c(1, 3, 6, 8, 9), ]
tempDat2

algW3 <- rbind(algW2, tempDat2)

with(algW3, table(st, YEAR))

algW4 <- droplevels(filter(algW3, st != "AQUE_6" &
                             st != "IVEE_3" &
                             st != "IVEE_4" &
                             st != "IVEE_5" &
                             st != "IVEE_6" &
                             st != "IVEE_7" &
                             st != "IVEE_8" ))
with(algW4, table(st, YEAR))

dat1 <- algW4
#################################
# Alpha - rich, div, even, abund
names(dat1)
dat1 <- divMetF(dat1, 5, 55)

# RENAME
datW <- dat1 %>% rename(site = SITE, year = YEAR)

qplot(dateR, rich, data = datW, 
      facets = ~site, geom = "point") + geom_smooth(method = "loess")

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
temp1 <- temp0[, 5:55] # MODIFY!
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
btsStudySub <- "SBCalg"
alphaSQM <- 40
gammaSQM <- NA
abundType <- "PercentCover"

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
  abundSQM = NA,
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
  abundSQM = NA,
  sqmOriginal = alphaSQM*alphDF$nReps[1], 
  nReps = NA
)

head(btsDFgamma)
btsDFgamma$subSiteID <- with(btsDFgamma, paste(studySub, site, sppCode, sep = "_"))

btsDF <- rbind(btsDFalpha, btsDFgamma)
head(btsDF)
tail(btsDF)

