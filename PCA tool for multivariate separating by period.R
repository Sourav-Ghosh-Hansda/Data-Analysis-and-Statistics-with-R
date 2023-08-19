if(!is.null(dev.list())) dev.off();rm(list = ls());cat("\014") # clear out the past 
library(dplyr)

setwd("C:/Users/dbrawn/Desktop/MA334 project")
Proj_data <-  read.csv("proportional_species_richness_Dan.csv") # you use V2 or V3
Proj_data$period <- as.factor(Proj_data$period) # must set categorical vars
Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)

# Proj_data_selected  <- Proj_data%>%filter(dominantLandClass=="7w") # Sea cliffs/hard coast, Wales
# Proj_data_selected  <- Proj_data%>%filter(dominantLandClass=="7e") # Sea cliffs/hard coast, England
# Proj_data_selected  <- Proj_data%>%filter(dominantLandClass=="6w") # Complex valley systems/table lands, Wales

#  Proj_data_selected <- Proj_data%>%filter(dominantLandClass=="8e") # Estuarine/soft coast/tidal rivers, England
# Proj_data_selected <- Proj_data%>%filter(dominantLandClass=="12e")  # Large river floodplains, flat plains, margins, E Anglia

# Proj_data_selected <- Proj_data%>%filter(dominantLandClass=="24s")  # Steep valley sides/intermediate mountain tops, W Highlands


BD_measures <- na.omit(Proj_data_selected) # need to remove NAs (or replace them) in a rotation
head(BD_measures[,2:12])
nrow(BD_measures)  # BD_measures still contains the periods in a single variable

# see ?prcomp the default here is the mean correct but not to scale 
pr.out=prcomp(BD_measures[,2:12]) # Principal Components 
pr.out$center  # gives the mean corrections; the "centers"
pr.out$scale  # not scaled
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2],col=BD_measures$period, cex=0.7,pch=16) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], BD_measures$period, cex=0.4, pos=4, col="red") # location labels

pr.out$rotation[,1:2] # the first two principal directions 

# watch out: the signs of the principal axes are not defined !
# to interpret the sign you could look at a dominating taxi group via boxplots:
plot(BD_measures$Carabids~BD_measures$period) # 7w & 7e :bad news for Isopods/Carabids (and us !)
plot(BD_measures$Isopods~BD_measures$period) # 8e & 12e :bad news for Isopods (and us !)
plot(BD_measures$Bees~BD_measures$period) # bees often increase but see 6w 
plot(BD_measures$Ladybirds~BD_measures$period) # RE 12e 