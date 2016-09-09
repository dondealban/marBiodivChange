#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Figure 2
# Author: Robin Elahi
# Date: 150629
#################################################

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# get data for analysis
source("./R/bts_dataPrep.R")

library(ggplot2)
library(nlme)
library(arm)
library(AICcmodavg)
library(plyr)

# plotting functions
source("./R/Standard Graphical Pars_RE.R")
source("./R/multiplotF.R")

# Random effects for LME models
source("./R/bts_lme_models.R")

# make sure the data are ordered by time within time-series (
# necessary for autocorrelation 
head(richDat$year0Z)
tail(divDat$year0Z)

############################################################
############################################################
###HMM WITH ESTIMATED SLOPES
###FOR RICHNESS AND DIVERSITY
############################################################
############################################################

# richness
richDat2 <- richDat[with(richDat, order(studyName, studySub, subSiteID, year0Z)), ]
# richDat2$Prediction <- relevel(richDat2$Prediction, ref = "none")
summary(richDat2$Prediction)

richGlobal <- lme(fixed = richLN ~ 1 +  
	year0Z*Scale + year0Z*Prediction +
	year0Z*Trophic + year0Z*initialRichLN + 
	year0Z*durationZ, 
	data = richDat2, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
plot(richGlobal)

# diversity
divDat2 <- divDat[with(divDat, order(studyName, studySub, subSiteID, year0Z)), ]
# reclassify 'neutral' prediction to 'none', because only one time-series with neutral
divDat2$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", divDat2$Prediction)))
# divDat2$PredictionOld <- relevel(divDat2$PredictionOld, ref = "none")
summary(divDat2$PredictionOld)

divGlobal <- lme(fixed = div ~ 1 +
	year0Z*Scale + year0Z*PredictionOld +
	year0Z*Trophic + year0Z*initialDiv + 
	year0Z*durationZ, 
	data = divDat2, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1(), weights = varExp())

# residual plots
plot(richGlobal, main = "ln(S)")
plot(divGlobal, main = "H'")

#######################################################
#######################################################
#Extract the fitted slopes for each time series to use for histograms

# get fitted values for each model
richDat2$fitted <- fitted(richGlobal)
divDat2$fitted <- fitted(divGlobal)

# use ddply to extract slope and se
slopeRich <- ddply(richDat2, .(subSiteID), summarize, 
			tempN = length(richLN), 
			meanS = mean(richLN),
			fitSlope = coefficients(lm(fitted ~ year0Z))[2]
			) 

slopeDiv <- ddply(divDat2, .(subSiteID), summarize, 
			tempN = length(div), 
			meanS = mean(div),
			fitSlope = coefficients(lm(fitted ~ year0Z))[2]
			) 

# merge slope data with original dataset by subsiteID
names(richDat2)
dim(richDat2)
dim(slopeRich)

tempRich <- merge(slopeRich, richDat2, by.x = "subSiteID")
tempDiv <- merge(slopeDiv, divDat2, by.x = "subSiteID")

# i just want one row for each subsiteID
slopeRich1 <- ddply(tempRich, .(subSiteID), head, n = 1)
slopeDiv1 <- ddply(tempDiv, .(subSiteID), head, n = 1)
# rename
slopeRich <- slopeRich1
slopeDiv <- slopeDiv1

# save the fitted slopes and export as csv
# (for Mark Vellend)
head(slopeRich)
richDF <- slopeRich %>% 
  dplyr::select(studyName, studySub, site, subSiteID, Scale, fitSlope)
head(richDF)
write.csv(richDF, "output/richness_slopes.csv")

divDF <- slopeDiv %>% 
  dplyr::select(studyName, studySub, site, subSiteID, Scale, fitSlope)
write.csv(divDF, "output/diversity_slopes.csv")

############################################################
############################################################
# 6 panel figure
with(slopeRich, table(Prediction))
richDrivers <- c(9, 22, 230, 41)

with(slopeDiv, table(PredictionOld))
divDrivers <- c(6, 157, 6)

ylab1 <- expression(paste("Change in richness [ln(S) ", yr^-1, "]"))
ylab2 <- expression(paste("Change in diversity [H' ", yr^-1, "]"))

xlab0 <- expression(paste("ln(S) ", yr^-1))
xlab1 <- expression(paste("H' ", yr^-1))

range(slopeRich$fitSlope)
range(slopeDiv$fitSlope)

set_graph_pars(ptype = "panel6")

hist(slopeRich$fitSlope, main = "", 
		xlab = ylab1, cex.axis = 1.2, cex.lab = 1.4, 
		col = "gray", breaks = 12, xlim = c(-0.06, 0.06), ylim = c(0,120))
abline(v = 0, lwd = 3, lty = 2, col = "black")
mtext(text = "A", side = 3, line = 1, adj = -0.1, cex = 1.4)
box()

plot(fitSlope ~ Prediction, data = slopeRich, 
			xlab = "Driver", ylab = ylab1, 
			cex.axis = 1.2, cex.lab = 1.4, 
			ylim = c(-0.06, 0.065))
abline(a = 0, b = 0, col = "darkgray", lty = 2, lwd = 2)
mtext(text = "C", side = 3, line = 1, adj = -0.1, cex = 1.4)
text(x = c(1:4), y = - 0.06, richDrivers, cex = 1.2)

plot(fitSlope ~ initialRichLN, data = slopeRich, 
			xlab = "Initial ln(S)", ylab = ylab1, 
			cex.axis = 1.2, cex.lab = 1.4,
			ylim = c(-0.06, 0.065))
abline(a = 0, b = 0, col = "darkgray", lty = 2, lwd = 2)
mtext(text = "E", side = 3, line = 1, adj = -0.1, cex = 1.4)

##2nd row
hist(slopeDiv$fitSlope, main = "", 
		xlab = ylab2, cex.axis = 1.2, cex.lab = 1.4, 
		col = "gray", breaks = 9, xlim = c(-0.06, 0.06), ylim = c(0,55))
abline(v = 0, lwd = 3, lty = 2, col = "black")
mtext(text = "B", side = 3, line = 1, adj = -0.1, cex = 1.4)
box()

plot(fitSlope ~ PredictionOld, data = slopeDiv, 
			xlab = "Driver", ylab = ylab2, 
			cex.axis = 1.2, cex.lab = 1.4, 
			ylim = c(-0.06, 0.06), names = c("negative", "none*", "positive"))
abline(a = 0, b = 0, col = "darkgray", lty = 2, lwd = 2)
mtext(text = "D", side = 3, line = 1, adj = -0.1, cex = 1.4)
text(x = c(1:3), y = - 0.06, divDrivers, cex = 1.2)

plot(fitSlope ~ initialDiv, data = slopeDiv, 
			xlab = "Initial H'", ylab = ylab2, 
			cex.axis = 1.2, cex.lab = 1.4,
			ylim = c(-0.06, 0.06))
abline(a = 0, b = 0, col = "darkgray", lty = 2, lwd = 2)
mtext(text = "F", side = 3, line = 1, adj = -0.1, cex = 1.4)

