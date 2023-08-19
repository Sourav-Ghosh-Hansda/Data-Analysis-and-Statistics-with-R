if(!is.null(dev.list())) dev.off()  # clear out the past 
rm(list = ls())
cat("\014")
library(dplyr)
library(ggplot2)
library(tidyr) # for spliting on the period see below
library(moments) # for calculating moments for skewness etc.
library(reshape2)
library(corrplot)
library(stargazer)
par(mfrow=c(1, 1)) 
#Load the data
Proj_data <-  read.csv("/Users/souravghoshhansda/Library/CloudStorage/OneDrive-UniversityofEssex/MA334/Project/folder/Materials/proportional_species_richness_V3.csv") # you use V2 or V3
Proj_data$period <- as.factor(Proj_data$period) # must set categorical vars
Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)
names(Proj_data)
#  count incidents (both periods) for each land classification for later selection
Proj_data%>%group_by(dominantLandClass)%>%count()%>%
  arrange(dominantLandClass)%>%print(n=45)
# select  7 randomly chosen predictors to form the trial eco_stat
all <- c(2:12)
eco_selected <- sample(all,size=7, replace = FALSE)
eco_selected <- c(2,3,4,5,6,10,12) #a particular troublesome case
eco_not_selected <- all[!(all%in%eco_selected)]
eco_names <- names(Proj_data[,2:12])
eco_selected_names <- names(Proj_data)[eco_selected]
eco_selected_names
# calculate the bio div measure over 7 taxinomic groups
mean_selected <- rowMeans(Proj_data[,eco_selected],na.rm=TRUE) # mean the 7 columns 
sum(is.na(mean_selected)) # check that there are no NAs in mean_selected
# add in the biodiversity measure which is the mean over 7 taxonomic groups
Proj_data_MA334 <- Proj_data%>%mutate(eco_status_7=mean_selected)
names(Proj_data_MA334)
# the data exploration phase 
table <- data.frame()
for(i in eco_selected){
  table <- rbind(table,
                 c(eco_names[i-1],
                   round(mean(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(median(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(min(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(max(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(sd(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(skewness(Proj_data_MA334[,i],na.rm = TRUE),digits = 2),
                   round(kurtosis(Proj_data_MA334[,i],na.rm = TRUE),digits = 2)
                 ))}
colnames(table) <- c("taxi_group","mean","median","min","max","sd","skewness","kurtosis")
table%>%arrange(mean,median,min,max,sd,skewness,kurtosis) # something more could be done here 
# extend data exploration; with correlations between continuous variables
names(Proj_data_MA334)
cont_vars <- Proj_data_MA334%>%select(c(eco_selected,13,14,16)) # includes   
names(cont_vars)
cormat <- round(x = cor(cont_vars,use="pairwise.complete.obs"), digits = 2)
# melt the correlation matrix
#melt(cormat)%>%mutate(R2 = value^2)%>%arrange(value)
melt(cormat)%>%mutate(R2 = value^2)%>%arrange(Var1,value)
#Plotting the relationships
plot(Proj_data_MA334$eco_status_7~Proj_data_MA334$ecologicalStatus)
cor(Proj_data_MA334$eco_status_7,Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Bees~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Bird~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Bryophytes~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Butterflies~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Carabids~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Macromoths~Proj_data_MA334$ecologicalStatus)
plot(Proj_data_MA334$Vascular_plants~Proj_data_MA334$ecologicalStatus)
# calculate the correlation matrix
corr_matrix <- cor(Proj_data_MA334[, eco_selected_names], use="complete.obs")
corrplot(corr_matrix, method = "color", type = "lower", tl.col = "black", tl.srt = 45, addCoef.col = "black")
#perform hypothesis test
names(Proj_data_MA334)
Proj_data_MA334_period <- Proj_data_MA334%>%select(Location,period,eco_status_7)
Proj_data_MA334_split <- Proj_data_MA334_period%>%pivot_wider(names_from =period,values_from=eco_status_7)
Proj_data_MA334_split <- Proj_data_MA334_split%>%mutate(BD7_change=Y00-Y70)
head(Proj_data_MA334_split)
hist(Proj_data_MA334_split$BD7_change)  # the distribution of the BD7 change 
BD7_change <- Proj_data_MA334_split%>%pull(BD7_change)
t.test(BD7_change,mu=0)  # t test with H0: mu=0
# comparing the two distributions of bio div based on 7 and 11 taxonomic groups 
par(mfrow=c(1, 1))  # divide graph area in 1 columns
qqplot(Proj_data_MA334$eco_status_7,Proj_data_MA334$ecologicalStatus)
abline(0,1,col="red")
# both cdfs together  and do a kolmogorov test H0: distributions are the same
# Calculate the empirical cumulative distribution functions (CDFs)
BD7_cdf <- ecdf(Proj_data_MA334$eco_status_7)
BD11_cdf <- ecdf(Proj_data_MA334$ecologicalStatus)
# Plot the CDFs
plot(BD11_cdf, col = "red", main = "Empirical CDFs of BD7 and BD11",
     xlab = "Proportional Species Richness", ylab = "CDF", xlim = c(0, 1), ylim = c(0, 1))
lines(BD7_cdf, col = "green")
# Add legends
legend("bottomright", legend = c("BD11", "BD7"), col = c("red", "green"), lty = 1)
ks.test(Proj_data_MA334$eco_status_7,Proj_data_MA334$ecologicalStatus)
# Simple linear regression 
# regressions of eco_status_7 against ecologicalstatus based on all 11
plot(Proj_data_MA334$eco_status_7~Proj_data_MA334$ecologicalStatus)
abline(0,1,col="red")
lin_mod <- lm(Proj_data_MA334$eco_status_7~Proj_data_MA334$ecologicalStatus)
lin_mod$coefficients
summary(lin_mod)
summary(lin_mod)$r.squared
abline(lin_mod,col="green")
legend("topleft", legend = c("Reference Line", "Fitted Regression Line"), col = c("red", "green"), lty = 1, cex = 0.8)
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")
qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")
# doing the same for each period report and differences 
Proj_data_MA334_Y70 <- Proj_data_MA334%>%filter(period=="Y70")
plot(Proj_data_MA334_Y70$eco_status_7~Proj_data_MA334_Y70$ecologicalStatus)
abline(0,1,col="red")
lm_mod_Y70 <- lm(Proj_data_MA334_Y70$eco_status_7~Proj_data_MA334_Y70$ecologicalStatus)
lm_mod_Y70$coefficients
summary(lm_mod_Y70)  
abline(lm_mod_Y70,col="green")
legend("topleft", legend = c("Reference Line", "Fitted Regression Line"), col = c("red", "green"), lty = 1, cex = 0.8)
# for later period 
Proj_data_MA334_Y00 <- Proj_data_MA334%>%filter(period=="Y00")
plot(Proj_data_MA334_Y00$eco_status_7~Proj_data_MA334_Y00$ecologicalStatus)
abline(0,1,col="red")
lm_mod_Y00 <- lm(Proj_data_MA334_Y00$eco_status_7~Proj_data_MA334_Y00$ecologicalStatus)
lm_mod_Y00$coefficients
summary(lm_mod_Y00) 
abline(lm_mod_Y00,col="green")
legend("topleft", legend = c("Reference Line", "Fitted Regression Line"), col = c("red", "green"), lty = 1, cex = 0.8)
# Multiple linear regression of BD4 on 7 taxinomical groups 
mean_selected <- rowMeans(Proj_data[,eco_not_selected ],na.rm=TRUE) # mean the rem 4 columns 
sum(is.na(mean_selected)) # check that there are no NAs in mean_selected
# add in the biodiversity measure which is the mean over 7 taxonomic groups
Proj_data_MA334 <- Proj_data_MA334%>%mutate(eco_status_4=mean_selected)
names(Proj_data_MA334)
# Create Training and Test data 
trainingRowIndex <- sample(1:nrow(Proj_data_MA334), 0.8*nrow(Proj_data_MA334))  # row indices for 80% training data
trainingData <- Proj_data_MA334[trainingRowIndex, ]  # model training data
testData  <- Proj_data_MA334[-trainingRowIndex, ]%>%na.omit # for test data remove NAs 
# Build the model on training data
lmMod_train <- lm(eco_status_4~.,
                  data=trainingData[c(eco_selected_names,"eco_status_4")],
                  na.action=na.omit,y=TRUE)
summary (lmMod_train)  # model summary
cor(lmMod_train$fitted.values,lmMod_train$y) # cor training data 
Eco_4_Pred <- predict(lmMod_train, testData) # predict to check model on test Data
cor(Eco_4_Pred,testData$eco_status_4)
plot(Eco_4_Pred~testData$eco_status_4)
abline(0,1,col="red")
# mis_fit_to_testData are the residuals for the train model fit to the test data 
mis_fit_to_testData <- testData$eco_status_4-Eco_4_Pred
plot(mis_fit_to_testData~Eco_4_Pred) # look for unwanted pattern in residuals
abline(0,0,col="red")
qqnorm(mis_fit_to_testData) # check for normality of residuals in prediction
qqline(mis_fit_to_testData,col="red")
backward <- step(lmMod_train, direction='backward', scope=formula(lmMod_train), trace=0)
backward$anova
backward$coefficients
forward <- step(lmMod_train, direction='forward', scope=formula(lmMod_train), trace=0)
forward$anova
forward$coefficients
#Analysing the backward-selected model's findings
summary(backward)
#Analysing the forward-selected model's findings
summary(forward)
#Open Analysis
# box plot comparisons for the two periods ignoring all other varaibles 
eco_status <- Proj_data_MA334%>%pull(eco_status_7)
eco_period <- Proj_data_MA334%>%pull(period)
plot(eco_status~eco_period)
# multiple linear regression BD7 against period, easting and northing 
mult_lin_mod <- lm(eco_status_7~.,
                   data=Proj_data_MA334[c("eco_status_7",
                                          "period","Easting","Northing")],
                   na.action = na.omit,y=TRUE)
summary(mult_lin_mod)
plot(mult_lin_mod$fitted.values~mult_lin_mod$y)
abline(0,1,col="red")
plot(jitter(fitted(mult_lin_mod)),residuals(mult_lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")
qqnorm(residuals(mult_lin_mod))
qqline(residuals(mult_lin_mod),col="red")
# compare the effect of each significant coefficient to that of period
mult_lin_mod$coefficients
as.numeric(mult_lin_mod$coefficients[3])*mean(Proj_data_MA334$Easting)
as.numeric(mult_lin_mod$coefficients[4])*mean(Proj_data_MA334$Northing)
# The following PCA method is an extension to the set book 
# PCA for visualizing the multi-dimensional spread of biodiversity values #######################
table(Proj_data_MA334_Y70$period); table(Proj_data_MA334_Y00$period) # check that these separate periods 
table(Proj_data_MA334_Y00$Location==Proj_data_MA334_Y70$Location) # check that Locations correspond between the two periods
eco_difference <- Proj_data_MA334_Y00[,eco_selected ]-Proj_data_MA334_Y70[,eco_selected ] # general differences between the two periods 
head(eco_difference)
# see ?prcomp the default here is the mean correct but not to scale 
pr.out=prcomp(na.omit(eco_difference)) # Principal Components 
pr.out$center  # gives the mean corrections the "centers"
pr.out$scale  # not scaled
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], Proj_data_MA334_Y00$dominantLandClass, cex=0.5, pos=4, col="red") # location labels
# label by location 
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], Proj_data_MA334_Y00$Location, cex=0.4, pos=4, col="red") # location labels
# label by eco increase 
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
BD_inc  <- Proj_data_MA334_Y00$eco_status_7-Proj_data_MA334_Y70$eco_status_7 # BD differences between the two periods 
text(pr.out$x[,1],pr.out$x[,2], round(BD_inc,2), cex=0.4, pos=4, col="red") # location labels
# label by a particular taxi group (if any) dominant in the first PC 
eco_selected <- c(2,3,4,5,6,10,12) # an example Bees !
eco_difference <- Proj_data_MA334_Y00[,eco_selected ]-Proj_data_MA334_Y70[,eco_selected ] # general differences between the two periods 
pr.out=prcomp(na.omit(eco_difference)) # Principal Components 
pr.out$rotation[,1:2] # print out first two principal axes
screeplot(pr.out, type="lines") # plot the variances in decreasing order
plot(pr.out$x[,1],pr.out$x[,2]) # scatter plot for first two principal components
text(pr.out$x[,1],pr.out$x[,2], round(eco_difference$Bees,2), cex=0.5, pos=4, col="red") # location labels