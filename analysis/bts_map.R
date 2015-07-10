#################################################
# Elahi et al. 2015. Recent trends in local-scale marine 
# biodiversity reflect community structure and human impacts
# Current Biology, in press June 2015

# Map for Figure 1
# Author: Robin Elahi
# Date: 150629
#################################################

library(rworldmap)
library(ggplot2)
library(gridExtra)
library(directlabels)
library(dplyr)

# rm(list=ls(all=TRUE)) 

siteList <- read.csv("./data/siteList_150617.csv")
names(siteList)
dim(siteList)
unique(siteList$site)
unique(siteList$studyName)

##############################
# drop the following studies
with(siteList, c(length(unique(studyName)), length(studyName)))

siteList <- siteList[siteList$studyName != "Keller" &
                       siteList$studyName != "Bebars" &
                       siteList$studyName != "Greenwood" &
                       siteList$studyName != "Sonnewald" &
                       siteList$studyName != "SwedFishTrawl" , ]
siteList <- droplevels(siteList)					
##############################
unique(siteList$site) # 191 levels
unique(siteList$studyName)

# create separate list of 2 studies with diversity only
siteList2 <- droplevels(filter(siteList, studyName == "Shin" |
                                 studyName == "Touzri"))
unique(siteList2$studyName)


# remove Shin, Touzri from first layer of points
siteList3 <- droplevels(filter(siteList, studyName != "Shin" &
                      studyName != "Touzri"))
unique(siteList3$studyName)
siteList <- siteList3

### For site list (55 studies)
dplyr::filter(siteList, studyName == "IMOS")
dplyr::filter(siteList, site == "TAR")

siteList$SiteIndex <- 1:length(siteList$site)
unique(siteList$SiteIndex)
unique(siteList$studyName)
siteList$Lat
siteList$Lat2 <- signif(siteList$Lat, 5)
siteList$Lat2

siteList$Long
siteList$Long2 <- signif(siteList$Long, 6)
siteList$Long2

### For site list2 (2 studies)
siteList2$SiteIndex <- 1:length(siteList2$site)
unique(siteList2$SiteIndex)
unique(siteList2$studyName)
siteList2$Lat
siteList2$Lat2 <- signif(siteList2$Lat, 5)
siteList2$Lat2

siteList2$Long
siteList2$Long2 <- signif(siteList2$Long, 6)
siteList2$Long2

# Get world map from rworldmap
worldmap <- map_data(map = "world")
# For interest you can take a look at the data for the worldmap
str(worldmap)

#######################
# Plot constrained world map with sites
latlimits <- c(-50, 75) 
longlimits <- c(-165, 170)  

no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())

sitemap <- ggplot(data = worldmap, aes(x = long, y = lat, group = group)) +
  geom_polygon(colour = "lightgray", fill = "lightgray") +
  coord_cartesian(xlim = longlimits, ylim = latlimits) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14)) + 
  no_grid

no_legend <- theme(legend.position = "none")

sitemap        

#######################
ULClabel <- theme(plot.title = element_text(hjust = -0.02, vjust = 2.1, size = rel(1.5)))

map1 <- sitemap + 
  geom_point(data = siteList, 
             aes(x = Long2, y = Lat2, group = studyName, color = studyName), 
             position = "jitter",  size = 3, alpha = 0.45) + 
  geom_point(data = siteList, 
             aes(x = Long2, y = Lat2, group = studyName, color = studyName), 
             position = "jitter",  size = 3, alpha = 1, shape = 21) + 
  geom_point(data = siteList2, 
             aes(x = Long2, y = Lat2, group = studyName, color = studyName), 
             position = "jitter",  size = 3, alpha = 0.45) + 
  geom_point(data = siteList2, 
             aes(x = Long2, y = Lat2, group = studyName, color = studyName), 
             position = "jitter",  size = 3, alpha = 1, shape = 21) + 
  no_legend + labs(title = "A") + ULClabel


