#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Funnel plots for Figure S2
# Author: Robin Elahi
# Date: 150629
#################################################

# get data for analysis
source("./R/bts_dataPrep.R")

# detach dplyr
detach("package:dplyr", unload = TRUE)

library(ggplot2)
library(nlme)
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
divDat2$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", divDat2$Prediction)))

divGlobal <- lme(fixed = div ~ 1 +
                   year0Z*Scale + year0Z*PredictionOld +
                   year0Z*Trophic + year0Z*initialDiv + 
                   year0Z*durationZ, 
                 data = divDat2, method = "REML", 
                 random =  list(rand1, rand2, rand4), 
                 correlation = corAR1(), weights = varExp())
plot(divGlobal)

# residual plots
plot(richGlobal, main = "ln(S)")
plot(divGlobal, main = "H'")

#######################################################
#######################################################
#Extract the fitted slopes for each time series to use for histograms
# get fitted values for each model

richDat2$fittedRich <- fitted(richGlobal)
divDat2$fittedDiv <- fitted(divGlobal)

# use ddply to extract slope and se
slopeRich <- ddply(richDat2, .(subSiteID, Prediction), summarize, 
                   tempN = length(richLN), 
                   meanS = mean(richLN),
                   fitSlope = coefficients(lm(fittedRich ~ year0Z))[2]
) 

unique(divDat2$PredictionOld)
slopeDiv <- ddply(divDat2, .(subSiteID, PredictionOld), summarize, 
                  tempN = length(div), 
                  meanS = mean(div),
                  fitSlope = coefficients(lm(fittedDiv ~ year0Z))[2]
) 

#################################
#################################
# FUNNEL PLOTS
#################################
#################################
# funnel plot for richness
lab0 <- expression(paste("Absolute value [ln(S) ", yr^-1, "]"))
lab1 <- expression(paste("Change in richness [ln(S) ", yr^-1, "]"))
lab2 <- expression(paste("Change in diversity [H' ", yr^-1, "]"))

# PUB QUALITY PLOTS
axesSpecs <- theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10))
no_grid <- theme(panel.grid.major = element_blank()) + 
theme(panel.grid.minor = element_blank())
plotSpecs <- theme_bw() + no_grid
regression <- geom_smooth(method = "lm", se = FALSE, size = .2)
yAxis <- scale_y_continuous(limits = c(0,1))
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 0.5, size = rel(1)))

slopeRich$Prediction
###

funnelRich <- ggplot(data = slopeRich, 
	aes(fitSlope, tempN, color = Prediction)) + 
	geom_point(size = 2, alpha = 0.1) +
	geom_point(size = 2, alpha = 1, pch = 21) +
	plotSpecs + axesSpecs + 
	geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
	ylab("Temporal replicates") + xlab(lab1) +
	scale_x_continuous(limits = c(-0.06, 0.06)) + 
	theme(legend.title = element_blank()) +
		scale_colour_manual(values = c("red", "darkgray", "black", "blue")) +
	labs(title = "C") + ULClabel + 
	theme(legend.position = "none")
funnelRich

levels(slopeRich$Prediction)
levels(slopeDiv$PredictionOld)

funnelDiv <- ggplot(data = slopeDiv, 
	aes(fitSlope, tempN, color = PredictionOld)) + 
	geom_point(size = 2, alpha = 0.1) +
	geom_point(size = 2, alpha = 1, pch = 21) +
	plotSpecs + axesSpecs + 
	geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
	ylab("Temporal replicates") + xlab(lab2) +
	scale_x_continuous(limits = c(-0.05, 0.05)) + 
	theme(legend.title = element_blank()) +
	scale_colour_manual(breaks = c("negative", "neutral", "positive"), 
		values = c("red", "black", "blue")) +
	labs(title = "D") + ULClabel + 
	theme(legend.justification = c(1, 1), 
		legend.position = c(0.4, 1)) +
	theme(legend.position = "none")
	
funnelDiv
		