#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Estimate (roughly) how many time series significantly
# positive or negative, in response to reviewer suggestion
# Author: Robin Elahi
# Date: 150629
#################################################
# Get fitted slopes from global model
# Get confidence interval for the year random effect
# Add the confidence intervals to fitted slopes, calculate % that don't overlap with zero

# rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# Source the slope data used for figure 2
source("./bts_Fig2.R")

head(slopeRich)
head(slopeDiv)

############################################################
############################################################
### RICHNESS CATERPILLAR PLOT
############################################################
############################################################
qplot(fitSlope, data = slopeRich)

# select model to use for error bars
mod <- richGlobal
summary(mod)
intervals(mod)$fixed
intervals(mod, which = 'var-cov')

# get estimate of the sd of year (random)
yearSD <- intervals(mod, which = 'var-cov')$reStruct$subSiteID["sd(year0Z", "est."]
# CI is 1.96*SD
yearCI <- 1.96*yearSD
yearCI # approximate 95% CI for the random effect

df1 <- slopeRich[with(slopeRich, order(fitSlope)), ]
df1$newRow <- seq(1:dim(df1)[1])

df1$betaUpper <- df1$fitSlope + yearCI
df1$betaLower <- df1$fitSlope - yearCI

sigUpper <- dim(df1[df1$betaLower > 0, ])[1]
sigLower <- dim(df1[df1$betaUpper < 0, ])[1]
sigUpper; sigLower
sigUpper/dim(df1)[1] # 16% sig pos
sigLower/dim(df1)[1] # 3% sig neg
sigPercentage <- (sigUpper + sigLower)/dim(df1)[1]; sigPercentage

# caterpillar plot
no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())
yAxis1 <- theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14))
plotSpecs <- theme_bw() + yAxis1 + xAxis1 + no_grid

ggplot(df1, aes(newRow, fitSlope)) + 
  geom_errorbar(aes(ymin = betaLower, ymax = betaUpper), 
                size = 0.03, width = 0.0, 
                color = "black") + geom_point(size = 1) + 
  geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") + 
  plotSpecs +
  xlab("Time-series") + ylab("Slope coefficient (year)") 

############################################################
############################################################
### DIVERSITY CATERPILLAR PLOT
############################################################
############################################################
qplot(fitSlope, data = slopeDiv)

# select model to use for error bars
mod <- divGlobal
summary(mod)
intervals(mod)$fixed
intervals(mod, which = 'var-cov')

# get estimate of the sd of year (random)
yearSD <- intervals(mod, which = 'var-cov')$reStruct$subSiteID["sd(year0Z", "est."]
# CI is 1.96*SD
yearCI <- 1.96*yearSD
yearCI # approximate 95% CI for the random effect

df1 <- slopeDiv[with(slopeDiv, order(fitSlope)), ]
df1$newRow <- seq(1:dim(df1)[1])

df1$betaUpper <- df1$fitSlope + yearCI
df1$betaLower <- df1$fitSlope - yearCI

sigUpper <- dim(df1[df1$betaLower > 0, ])[1]
sigLower <- dim(df1[df1$betaUpper < 0, ])[1]
sigUpper; sigLower
sigUpper/dim(df1)[1] # 6% sig pos
sigLower/dim(df1)[1] # 3% sig neg
sigPercentage <- (sigUpper + sigLower)/dim(df1)[1]; sigPercentage

# caterpillar plot
ggplot(df1, aes(newRow, fitSlope)) + 
  geom_errorbar(aes(ymin = betaLower, ymax = betaUpper), 
                size = 0.03, width = 0.0, 
                color = "black") + geom_point(size = 1) + 
  geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") + 
  plotSpecs +
  xlab("Time-series") + ylab("Slope coefficient (year)") 

