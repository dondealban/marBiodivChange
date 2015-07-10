#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Figure S2
# Author: Robin Elahi
# Date: 150629
#################################################

library(ggplot2)

source("./R/multiplotF.R")

# load the simulated data
pValsLongRich <- read.csv("./output/pValsLong10reps_150620keep.csv", header=TRUE, 
                          na.strings="NA")

pValsLongDiv <- read.csv("./output/pValsLongDiversity_150620keep.csv", header=TRUE, 
                         na.strings="NA")

dim(pValsLongRich)
# 36 classifications x 10 reps each x 1 row per 4 p-values (neg, pos, neut, overall)
36 * 10 * 4 # but only 1400 rows - missing '65' (x-axis)

dim(pValsLongDiv)
13 * 10 * 3 # appropriate number

pValsLongRich65 <- read.csv("./output/pValsLong10reps_65_150620keep.csv", header=TRUE, 
                          na.strings="NA")
dim(pValsLongRich65)
pValsLongRich2 <- rbind(pValsLongRich, pValsLongRich65)
unique(pValsLongRich2$predictionChanges)

# PUB QUALITY PLOTS
axesSpecs <- theme(axis.title = element_text(size = 12), 
                   axis.text = element_text(size = 10))
no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())
plotSpecs <- theme_bw() + no_grid
regression <- geom_smooth(method = "lm", se = FALSE, size = .2)
yAxis <- scale_y_continuous(limits = c(0,1))
ULClabel <- theme(plot.title = element_text(hjust = -0.07, vjust = 0.05, size = rel(1)))


# without facets, manual colors
divPlot <- ggplot(pValsLongDiv, aes(nChanges, value, 
                                    color = variable)) + 
  geom_point(size = 2, alpha = 0.1) + 
  geom_point(size = 2, alpha = 0.8, pch = 21) +
  plotSpecs + axesSpecs + 
  ylab("P-value") + xlab("Number of alternative classifications") + 
  geom_smooth(se = FALSE, size = 1) + 
  scale_x_continuous(breaks = seq(0,70, by = 5)) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(breaks = c("positive", "negative", "overall"), 
                      values = c("red", "black", "blue")) +
  labs(title = "B") + ULClabel +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  geom_hline(yintercept = 0.05, color = "green", linetype = "dashed", size = 1.5)

divPlot
7/13 # 54%
############################################################
############################################################
# plot richness version
richPlot <- ggplot(pValsLongRich2, aes(nChanges, value, 
                                      color = variable)) + 
  geom_point(size = 2, alpha = 0.1) + 
  geom_point(size = 2, alpha = 0.8, pch = 21) +
  plotSpecs + axesSpecs + 
  ylab("P-value") + xlab("Number of alternative classifications") + 
  geom_smooth(se = FALSE, size = 1) + 
  scale_x_continuous(breaks = seq(0,70, by = 5)) +
  theme(legend.title = element_blank()) +
  labs(title = "A") + ULClabel +
  scale_colour_manual(breaks = c("positive", "neutral", "negative", "overall"), 
                      values = c("red", "darkgray", "black", "blue")) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  geom_hline(yintercept = 0.05, color = "green", linetype = "dashed", size = 1.5)

richPlot
23/72 # 32% 
############################################################
############################################################
# source funnel plots
source("./bts_funnelPlot.R")

############################################################
############################################################
# FINAL FIGURE S2

pdf("./cb_FigS2.pdf", 6.9, 6.9)
multiplot(richPlot, funnelRich, divPlot, funnelDiv, cols = 2)
dev.off()	


