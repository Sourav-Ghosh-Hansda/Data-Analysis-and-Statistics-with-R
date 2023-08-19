if(!is.null(dev.list())) dev.off()  # clear out the past 
rm(list = ls())
cat("\014")

choose(11,7) # the number of separate of combinations of 7 from 11 taxonomic groups
library(dplyr)

setwd("/Users/souravghoshhansda/Library/CloudStorage/OneDrive-UniversityofEssex/MA334/Project/folder/Materials")
Proj_data <-  read.csv("proportional_species_richness.csv")
Proj_data$period <- as.factor(Proj_data$period) # must set categorical vars
Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)

names(Proj_data)
str(Proj_data)

# select  randomly chosen predictors to form the trial eco_stat

taxi_n <- 7
eco_selected <- sample(c(2:12),size=taxi_n, replace = FALSE)
eco_selected
eco_selected_names <- names(Proj_data)[eco_selected]
names(Proj_data[,eco_selected])

# selecting just the chosen 7 and the dom land class, eco measure and period
Var_demo <- Proj_data%>%select(all_of(eco_selected),dominantLandClass,
                               ecologicalStatus,period)
head(Var_demo)

mean_selected <- rowSums(Var_demo[1:7],na.rm=TRUE)/7 # mean the 7 columns 
sum(Var_demo[472,1:7])/7; mean_selected[472]# confirms the row  sums in a few cases as a check

par(mfrow=c(1, 2))  # divide graph area in 2 columns
hist(mean_selected);hist(Var_demo$ecologicalStatus)

par(mfrow=c(1, 1)) 
names(Var_demo)
hist(Var_demo$Isopods)
summary(Var_demo$Isopods)
mean(Var_demo$Isopods,trim=0.1)

# some examples of the use of the dplyr package
Var_demo%>%select(Isopods)%>%summary()

eco_status <- Var_demo%>%pull(ecologicalStatus)
eco_period <- Var_demo%>%pull(period)
boxplot(eco_status)
plot(eco_status~eco_period)
summary(lm(eco_status~eco_period)) # all to come on MA334 !!!
# compare to the difference in the means for the two periods
mean(Var_demo%>%filter(period=="Y70")%>%pull(ecologicalStatus))
mean(Var_demo%>%filter(period=="Y00")%>%pull(ecologicalStatus))

#lets look at the worse and best eco wise for both periods combined
hist(Var_demo$ecologicalStatus)
Qs <- quantile(Var_demo$ecologicalStatus,probs=seq(0,1,0.1))
str(Qs)

Best_eco <- Var_demo%>%filter(ecologicalStatus>Qs[10])%>%
  group_by(dominantLandClass)%>%count()%>%arrange(desc(n))

Worst_eco <- Var_demo%>%filter(ecologicalStatus<Qs[2])%>% 
  group_by(dominantLandClass)%>%count()%>%arrange(desc(n))


par(mfrow=c(1, 1))  # divide graph area in 2 columns
extreme_eco <- Var_demo%>%
  filter(ecologicalStatus>Qs[10])%>%
  select(ecologicalStatus)%>%pull(ecologicalStatus)
extreme_LC <- Var_demo%>%
  filter(ecologicalStatus>Qs[10])%>%
  select(period)%>%pull(period)
boxplot(extreme_eco)
plot(extreme_eco~extreme_LC)

extreme_eco <- Var_demo%>%
  filter(ecologicalStatus<Qs[2])%>%
  select(ecologicalStatus)%>%pull(ecologicalStatus)
extreme_LC <- Var_demo%>%
  filter(ecologicalStatus<Qs[2])%>%
  select(period)%>%pull(period)
boxplot(extreme_eco)
plot(extreme_eco~extreme_LC)

# there follows some analysis by land classification and period

mean_LC <- Var_demo%>%group_by(dominantLandClass)%>%
  summarise(LC_mean =mean(ecologicalStatus))
hist(mean_LC$LC_mean)

Var_demo%>%group_by(dominantLandClass,period)%>%
           count()%>%arrange(desc(n))

Var_demo%>%group_by(dominantLandClass,period)%>%
  summarise(LC_mean =mean(ecologicalStatus))%>%print(n=90)

library(tidyr) # for spliting on the period 
 eco_changes <- Var_demo%>%group_by(dominantLandClass,period)%>%
  summarise(LC_mean =mean(ecologicalStatus))%>%
pivot_wider(names_from = period, values_from = LC_mean, values_fill = 0)

 Var_demo%>%group_by(dominantLandClass,period)%>%
   summarise(LC_mean =mean(ecologicalStatus))%>%
   pivot_wider(names_from = period, values_from = LC_mean, values_fill = 0)%>%
   mutate(eco_change=Y00-Y70)
 
 Var_demo%>%group_by(dominantLandClass,period)%>%
   summarise(LC_mean =mean(ecologicalStatus),.groups = 'drop')%>%
   pivot_wider(names_from = period, values_from = LC_mean, values_fill = 0)%>%
   mutate(eco_change=Y00-Y70)%>%
   arrange(desc(eco_change))%>%print(n=45)  # note 33 out of 45 eco_change< 0
 
 