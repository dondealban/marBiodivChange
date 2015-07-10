#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Figure S1
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

# plotting functions
source("./R/Standard Graphical Pars_RE.R")
source("./R/multiplotF.R")

# Random effects for LME models
source("./R/bts_lme_models.R")

# make sure the data are ordered by time within time-series (
# necessary for autocorrelation 
head(divDat$year0Z)
tail(divDatSub$year0Z)


############################################################
############################################################
###FIGURE S1A: H' V YEAR
############################################################
############################################################
# NLME

ggDat <- divDat

simple <- lme(fixed = div ~ I(year0Z + 1990), 
	data = ggDat, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1(), weights = varExp())

intervals(simple)
summary(simple)

int <- fixed.effects(simple)[1]
slope <- fixed.effects(simple)[2]
int;slope
dim(ggDat)

ggDat$simple.fitted <- predict(simple)
ggDat$simple.resids <- with(ggDat, simple.fitted - div)
head(ggDat)
par(mfrow = c(1,2))
plot(simple.fitted ~ year0Z, ggDat)
plot(simple.resids ~ year0Z, ggDat)
year <- ggDat$year0 + 1962

# ggplot details
no_legend <- theme(legend.position = "none")

theme1 <- theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 14)) 

ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 0, size = rel(1.5)))

no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())

# Points colored by studyName, but fitted lines by subSiteID
# With abundance, solid lines; w/o = dashed lines

divP <- ggplot(ggDat, aes(year0 + 1962, div)) + theme1 + no_legend + 
  geom_point(aes(color = studyName), size = 1, alpha = 0.8) + 
  geom_smooth(data = subset(ggDat, AbundYes == 1), 
              method = "lm", se = F, size = 0.15, 
              aes(group = subSiteID, color = studyName)) + 
  geom_smooth(data = subset(ggDat, AbundYes == 0), 
              method = "lm", se = F, size = 0.15, linetype = "dashed", 
              aes(group = subSiteID, color = studyName)) 

FigS1A <- divP + xlab("Year") + ylab("Diversity (H')") +
  labs(title = "A") + ULClabel + 
  geom_abline(intercept = int, slope = slope, color = "black", size = 0.3) + no_grid

FigS1A

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
###TESTING HYPOTHESES FROM TABLE 1
############################################################
############################################################
lmeDat <- divDat
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")
# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")

# Use Gelman's to better match the other inputs

# Rescale duration to visualize coefficient better
#lmeDat$RS.durationZ <- scale(lmeDat$duration)[,1]
lmeDat$RS.durationZ <- rescale(lmeDat$duration)

# Rescale year to match duration
# lmeDat$RS.year0Z <- scale(lmeDat$year0)[,1]
lmeDat$RS.year0Z <- rescale(lmeDat$year0)

# How do the standardized inputs compare?
lmeDat %>% summarise(divSD = sd(div), 
                     initialDivSD = sd(initialDiv), 
                     yearSD = sd(RS.year0Z), 
                     durationSD = sd(RS.durationZ), 
                     driverSD = sd(PredictionOld),
                     trophicSD = sd(Trophic), 
                     scaleSD = sd(Scale))

cm1S <- lme(fixed = div ~ 1 +  
              RS.year0Z*Scale + RS.year0Z*PredictionOld +
              RS.year0Z*Trophic + RS.year0Z*initialDiv + 
              RS.year0Z*RS.durationZ, 
            data = lmeDat, method = "REML", 
            random =  list(rand1, rand2, rand4), 
            correlation = corAR1(), weights = varExp())

summary(cm1S)$tTable

cm4S <- update(cm1S, ~ . - RS.year0Z:Trophic)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)

bestModSet <- list(cm1S, cm4S, cm6S)
fullModel <- cm1S
bestNames <- paste("bestMod", c(1, 4, 6), sep = " ")

#######################################################
# write a loop to extract the rowname, coefficient, and upper and lower CL
# save the rownames
coefNames <- names(fixed.effects(fullModel))
coefNames
length(coefNames)
AllReps <- coefNames[10:16] # include only interactions
AllReps #
N <- length(AllReps); N
yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, 
	parm = "RS.year0Z", exclude = list(AllReps), second.ord = FALSE)

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
yearZname <- "year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, 
                    second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative","Year * Driver-Positive", "Year * Trophic Level 2", "Year * Trophic Level 3", "Year * Initial H'", "Year * Duration")
loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)
loopDat4$ci <- as.numeric(loopDat4$ci)

fullBetas <- loopDat4
divBetas <- fullBetas

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
###including abundance
############################################################
############################################################
lmeDat <- divDatSub
# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")

# Rescale duration to visualize coefficient better
lmeDat$RS.durationZ <- rescale(lmeDat$duration)

# Rescale year to match duration
lmeDat$RS.year0Z <- rescale(lmeDat$year0)

cm1S <- lme(fixed = div ~ 1 +  
              RS.year0Z*Scale + RS.year0Z*PredictionOld +
              RS.year0Z*Trophic + RS.year0Z*initialDiv + 
              RS.year0Z*RS.durationZ + abundCS, 
            data = lmeDat, method = "REML", 
            random =  list(rand1, rand2, rand4), 
            correlation = corAR1(), weights = varExp())

cm2S <- update(cm1S, ~ . - RS.year0Z:Scale)
cm5S <- update(cm1S, ~ . - RS.year0Z:initialDiv)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)

bestModSet <- list(cm1S, cm2S, cm5S, cm6S)
fullModel <- cm1S
bestNames <- paste("bestMod", c(1, 2, 5, 6), sep = " ")

#######################################################
coefNames <- names(fixed.effects(fullModel))
coefNames
length(coefNames)
AllReps <- coefNames[11:17] # include only interactions
AllReps #
N <- length(AllReps); N

yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, parm = "RS.year0Z", 
                    exclude = list(AllReps), second.ord = FALSE)

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
yearZname <- "year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, 
                    second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative","Year * Driver-Positive", 
             "Year * Trophic Level 2", "Year * Trophic Level 3", "Year * Initial H'", 
             "Year * Duration")
loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)
loopDat4$ci <- as.numeric(loopDat4$ci)

divBetasAbund <- loopDat4

############################################################
############################################################
###HIERARCHICAL MIXED EFFECTS MODELS: DIVERSITY
### REDUCED DATASET, WITHOUT ABUNDANCE IN MODELS
############################################################
############################################################
lmeDat <- divDatSub
# only one neutral prediction - lme won't run - need to use old classification
lmeDat$PredictionOld <- droplevels(as.factor(gsub("neutral", "none", lmeDat$Prediction)))
lmeDat$PredictionOld <- relevel(lmeDat$PredictionOld, ref = "none")

# Rescale duration to visualize coefficient better
lmeDat$RS.durationZ <- rescale(lmeDat$duration)

# Rescale year to match duration
lmeDat$RS.year0Z <- rescale(lmeDat$year0)

cm1S <- lme(fixed = div ~ 1 +  
              RS.year0Z*Scale + RS.year0Z*PredictionOld +
              RS.year0Z*Trophic + RS.year0Z*initialDiv + 
              RS.year0Z*RS.durationZ, 
            data = lmeDat, method = "REML", 
            random =  list(rand1, rand2, rand4), 
            correlation = corAR1(), weights = varExp())
summary(cm1S)$tTable

cm2S <- update(cm1S, ~ . - RS.year0Z:Scale)
cm5S <- update(cm1S, ~ . - RS.year0Z:initialDiv)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)

bestModSet <- list(cm1S, cm2S, cm5S, cm6S)
fullModel <- cm1S
bestNames <- paste("bestMod", c(1, 2, 5, 6), sep = " ")

#######################################################
# write a loop to extract the rowname, coefficient, and upper and lower CL
# save the rownames
coefNames <- names(fixed.effects(fullModel))
coefNames
length(coefNames)
AllReps <- coefNames[10:16] # include only interactions
AllReps #
N <- length(AllReps); N

yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, parm = "RS.year0Z", 
                    exclude = list(AllReps), second.ord = FALSE)

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
yearZname <- "year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative","Year * Driver-Positive", 
             "Year * Trophic Level 2", "Year * Trophic Level 3", "Year * Initial H'", 
             "Year * Duration")
loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)
loopDat4$ci <- as.numeric(loopDat4$ci)

divBetasReduced <- loopDat4

############################################################
############################################################
###MULTI-PLOT FIGURE
###DIV V YEAR AND MODEL AVERAGED COEFFICIENTS
############################################################
divBetas
divBetasReduced
divBetasAbund

# use these values to scale the size of the points
divScaleFull <- dim(divDat)[1]
divScaleRed <- dim(divDatSub)[1]

Betas <- rbind(divBetas, divBetasReduced, divBetasAbund)

Betas$Dataset <- c(rep("Full", 8), rep("Reduced", 16))
Betas$Models <- c(rep("No abundance", 16), rep("With abundance", 8))
Betas$DM <- c(rep("Full", 8), rep("Reduced", 8), rep("Reduced, abundance predictor", 8))

# order the values
divBetas[order(divBetas$modAvgBeta), ]
unique(Betas$xlabels)

# this will reorder based on fullBetas above
Betas$xlab2 <- factor(Betas$xlabels, levels = rev(c("Year * Initial H'", 
                                                    "Year * Duration",
                                                    "Year * Scale-Gamma", 
                                                    "Year * Driver-Negative",
                                                    "Year * Trophic Level 2",  
                                                    "Year * Trophic Level 3",
                                                    "Year * Driver-Positive", 
                                                    "Year")))

Betas
Betas <- droplevels(Betas)

# plot
no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())

yAxis1 <- theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14))
dodge <- position_dodge(width = -0.7)
plotSpecs <- theme1 + yAxis1 + xAxis1 

divScaleFull;divScaleRed

allBetas <- ggplot(Betas, aes(x = xlab2, y = modAvgBeta, group = DM, color = Models, size = Dataset)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") + plotSpecs + 
  geom_point(shape = 19, position = dodge) + 
  scale_size_manual(values = c(2.906*1.2, 2.586*1)) +
  geom_errorbar(width = 0, aes(ymin = lowerCL, ymax = upperCL), position = dodge, size = 0.8) + 
  scale_color_manual(values = c("black", "darkgray")) +
  coord_flip() + ylab("Standardized coefficients") +
  no_grid

ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 0, size = rel(1.5)))

mocFigDiv <- allBetas + theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  labs(title = "B") + ULClabel 
# mocFigDiv

############################################################
############################################################
### FIGURE S1, 2 panels
############################################################
############################################################

multiplot(FigS1A, mocFigDiv, cols = 2)
