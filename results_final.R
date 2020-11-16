# Script for presenting my results 

library(tidyverse)
library(ggplot2)  
library(grid)
library(viridis) 
library(raster)
library(rnaturalearth)
library(lme4)

worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')


# wm <- borders("world", colour = "gray50", fill = "gray50") 

# Part for inloading of data

analysis_stats <- readRDS("data//analysis//analysis_stats.rds")

# Filtering out models with AUC below 0.75

analysis_stats <- analysis_stats  %>%
  filter((AUC_mean) > 0.75)




# In case some ranges have completely contracted they will have NaN values for elevations, centroids, and so on.
# I overwrite Na- values with the values from the preceding time period
# overwriting na-values in the 2070 columns with data from 2050 

for (i in 1:nrow(analysis_stats)){
  if(is.na(analysis_stats$Centroid_x_2070[i]) == TRUE){
    analysis_stats$Centroid_x_2070[i] <- analysis_stats$Centroid_x_2050[i]
  }
  if(is.na(analysis_stats$Centroid_y_2070[i]) == TRUE){
    analysis_stats$Centroid_y_2070[i] <- analysis_stats$Centroid_y_2050[i]
  }
  if(is.na(analysis_stats$Mean_Elevation_2070[i]) == TRUE){
    analysis_stats$Mean_Elevation_2070[i] <- analysis_stats$Mean_Elevation_2050[i]
  }
  
}


dat <- read.csv("data//Plots//dat_12_10_thinned.csv", header = TRUE)
dat <- dat %>%
  filter(!is.na(elevationData))

ggplot(data = dat, aes( x = decimallongitude, y = decimallatitude)) +
  worldmap +
  geom_point()



ggplot(data = dat) +
  coord_fixed(ratio = 1) +
  geom_sf(data = worldmap) +
  geom_point(aes(x = decimallongitude, y = decimallatitude, color = (elevationData/1000))) +
  # geom_text(data = long_data, aes(x = Elevation, y = 5, label = Elevation)) +
  # geom_hline(yintercept = 35) +
  scale_color_viridis(discrete = FALSE, name = "Elevation (km)") +
  theme(legend.position="bottom") +
  
  # xlim(72,143) +
  # ylim(-70,70) +
  xlab("Longitude") +
  ylab("Latitude")





#### Making sf dataframe with coordnates for centroids

centroids <- analysis_stats[5:10]
dists_2050 <- pointDistance(centroids[1:2], centroids[3:4], lonlat = TRUE) / 1000
dists_2070 <- pointDistance(centroids[3:4], centroids[5:6], lonlat = TRUE) / 1000
# distances for doubling back: centroids moving back towards their present location
dists_db <- pointDistance(centroids[1:2], centroids[5:6], lonlat = TRUE) / 1000
total_distances <- dists_2050 + dists_2070

## looking at doubling back
db_values <- dists_db / total_distances


## calculating polewards shift

# delta_N_2050 <-  abs(analysis_stats$Centroid_y_2050 - analysis_stats$Centroid_y_Present)
# delta_N_2070 <-  abs(analysis_stats$Centroid_y_2070 - analysis_stats$Centroid_y_2050)
# delta_N_total <- delta_N_2050 + delta_N_2070  
#   
delta_N_2050 <- ifelse(analysis_stats$Centroid_y_Present > 0, analysis_stats$Centroid_y_2050 - analysis_stats$Centroid_y_Present, 
                       - (analysis_stats$Centroid_y_2050 - analysis_stats$Centroid_y_Present)
)

delta_N_2070 <- ifelse(analysis_stats$Centroid_y_2050 > 0, analysis_stats$Centroid_y_2070 - analysis_stats$Centroid_y_2050, 
                       - (analysis_stats$Centroid_y_2070 - analysis_stats$Centroid_y_2050)
)

delta_N_total <- delta_N_2050 + delta_N_2070  

dist_N_2050 <- pointDistance(centroids[1:2], cbind(centroids[c(1,4)]), lonlat = TRUE) / 1000
dist_N_2070 <- pointDistance(cbind(centroids[c(1,4)]), cbind(centroids[c(1,6)]), lonlat = TRUE) / 1000

analysis_stats <- cbind(analysis_stats, dists_2050, dists_2070, db_values, delta_N_2050, delta_N_2070, delta_N_total, dist_N_2050,dist_N_2070)

analysis_stats$Elevation_Group <- as.factor(ifelse(analysis_stats$Mean_Elevation_2070 > 3000, 1,0))



## calculating extinction risk as described in the PNAS article

extinction_2050 <- 1- (analysis_stats$Area_2050/analysis_stats$Present_Area)^0.25
extinction_2050 <- ifelse(extinction_2050 < 0, 0, extinction_2050)
extinction_2050 <- round(extinction_2050, 2)

extinction_2050_overlap <- 1- (analysis_stats$Area_overlap_2050/analysis_stats$Present_Area)^0.25
extinction_2050 <- ifelse(extinction_2050 < 0, 0, extinction_2050)
extinction_2050 <- round(extinction_2050, 2)


extinction_2070 <- 1- (analysis_stats$Area_2070/analysis_stats$Present_Area)^0.25
extinction_2070 <- ifelse(extinction_2070 < 0, 0, extinction_2070)
extinction_2070 <- round(extinction_2070, 2)

extinction_2070 <- 1- (analysis_stats$Area_2070/analysis_stats$Present_Area)^0.25
extinction_2070 <- ifelse(extinction_2070 < 0, 0, extinction_2070)
extinction_2070 <- round(extinction_2070, 2)




###  
long_data <- rbind(
  transmute(analysis_stats, 
            Period = as.factor("Present"),
            Species,
            Area = Present_Area,
            Centroid_X = Centroid_x_Present,
            Centroid_Y = Centroid_y_Present,
            Elevation = Present_Mean_Elevation,
            Contraction = 0,
            Distance = 0,
            Overlap = 1,
            Shift_Accelerate = as.factor(ifelse(dists_2050 - 0 > 0, 1,-1)),
            Double_Back = 0,
            Double_Back_Relative = 0,
            PoleShift = 0,
            Elevation_Group,
            PoleShift_km = 0,
            Extinction = 0
  ),
  
  transmute(analysis_stats, 
            Period = as.factor("2050"),
            Species,
            Area = Area_2050,
            Centroid_X = Centroid_x_2050,
            Centroid_Y = Centroid_y_2050,
            Elevation = Mean_Elevation_2050,
            Contraction = Range_Contraction_2050,
            Distance = dists_2050,
            Overlap = Area_overlap_2050 / Present_Area,
            Shift_Accelerate = as.factor(ifelse(dists_2070 - dists_2050 > 0, 1,-1)),
            Double_Back = 0,
            Double_Back_Relative = 0,
            PoleShift = delta_N_2050,
            Elevation_Group,
            PoleShift_km = dist_N_2050,
            Extinction = extinction_2050
  ),
  
  transmute(analysis_stats, 
            Period = as.factor("2070"),
            Species,
            Area = Area_2070,
            Centroid_X = Centroid_x_2070,
            Centroid_Y = Centroid_y_2070,
            Elevation = Mean_Elevation_2070,
            Contraction = Range_Contraction_2070,
            Distance = dists_2070,
            Overlap = Area_overlap_2070 / Present_Area,
            Shift_Accelerate = as.factor(0),
            Double_Back = dists_db,
            Double_Back_Relative = db_values,
            PoleShift = delta_N_2070,
            Elevation_Group,
            PoleShift_km = dist_N_2070,
            Extinction = extinction_2070
  )
  
) 

long_data$Extinction_Overlap <- round(1-(long_data$Overlap)^0.25,2)

## saving long_data

# saveRDS(long_data, file = "data//analysis//analysis_data.rds")
# write.csv(long_data, file = "data//analysis//analysis_data.csv")


##### Statistics stuff

# lm_2050 <- lm(dists_2050 ~ Present_Mean_Elevation + Present_Area + Centroid_y_Present + Centroid_y_Present*Present_Mean_Elevation +  
#                 Present_Mean_Elevation*Present_Area, data = analysis_stats)
# lm_2070 <- lm((dists_2070) ~ Mean_Elevation_2050 + Area_2050 + Centroid_y_2050 + Centroid_y_2050*Mean_Elevation_2050 + 
#                 Mean_Elevation_2050*Area_2050, data = analysis_stats)
# 
# 
# lm_2050 <- lm(Distance ~ Elevation + Area + Centroid_Y + Elevation*Centroid_Y, data = long_data[long_data$Period != "Present",])
# 
# lm_2050$coefficients
# lm_2070$coefficients
# 
# anova(lm_2050)
# summary(lm_2050)
# 
# anova(lm_2070)
# summary(lm_2070)
# 
# 

# trying out some glmer stuff

glm_test <- glm(Elevation_Group ~ Area + Centroid_Y, 
                data = long_data[long_data$Period == 2050,],
                family = "binomial")





#################### 
#
#  RANGE SIZE CHANGE PLOTS
#
##################  

# histogram of area changes for 2050 and 2070  
ggplot(data = long_data[long_data$Period != "Present",], aes(x = Contraction, color = Period)) +
  geom_histogram(aes(fill = Period), binwidth = 0.1, position = "dodge", color = "black") +
  xlab("Area change, km2")+
  ylab("Number of species")



## Scatter plot for area in relationship with future area changes
ggplot(data = long_data, aes(y = Contraction, x = log(Area+1), color = Elevation_Group, shape = Period)) +
  # coord_fixed(ratio = 1) +
  geom_point(alpha = 1, size = 2) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = Species)) +
  # theme(legend.position="none") +
  # scale_color_viridis(discrete = TRUE, option = "D") +
  # theme_dark() + 
  # ggtitle("Area changes") +
  scale_color_discrete(name = "Elevation Group", labels = c("Lowland", "Highland")) +
  # guides(color=guide_legend(title="Elevation Group")) +
  xlab("log(Present Area + 2) km^2") +
  ylab("Future Area Change")

# Scatter plot for contractions (x-axis) in relation with elevation (y-axis) for 2050 and 2070.  
ggplot(data = long_data[long_data$Period != "Present",] , aes( x = Contraction, y = Elevation)) +
  geom_point(aes(color = Period)) +
  # geom_line(aes(group = Species)) +
  xlab("Area change") +
  ylab("Elevation (m)")

# statistical test for the figure. see if contractions are different for the two elevation groups at both periods.
# use wilcox test. 

# first check if there is a contraction difference between periods  
wilcox.test((long_data$Contraction[long_data$Period == 2050]),
            (long_data$Contraction[long_data$Period == 2070]), 
            paired = FALSE)

# data:  (long_data$Contraction[long_data$Period == 2050]) and (long_data$Contraction[long_data$Period == 2070])
# W = 8012, p-value = 0.1313
# alternative hypothesis: true location shift is not equal to 0


wilcox.test((long_data$Contraction[long_data$Period == 2050 & long_data$Elevation_Group == "0"]),
            (long_data$Contraction[long_data$Period == 2050 & long_data$Elevation_Group == "1"]), 
            paired = FALSE)

# data:  (long_data$Contraction[long_data$Period == 2050 & long_data$Elevation_Group == "0"]) and (long_data$Contraction[long_data$Period == 2050 & long_data$Elevation_Group == "1"])
# W = 703, p-value = 0.1954
# alternative hypothesis: true location shift is not equal to 0

wilcox.test((long_data$Contraction[long_data$Period == 2070 & long_data$Elevation_Group == "0"]),
            (long_data$Contraction[long_data$Period == 2070 & long_data$Elevation_Group == "1"]), 
            paired = FALSE, alternative = "two.sided")
# 
# data:  (long_data$Contraction[long_data$Period == 2070 & long_data$Elevation_Group == "0"]) and (long_data$Contraction[long_data$Period == 2070 & long_data$Elevation_Group == "1"])
# W = 868, p-value = 0.958
# alternative hypothesis: true location shift is not equal to 0


## Boxplot for area changes in 2050 and 2070.  
ggplot(data = long_data, aes(x=Period, y = Contraction, fill = Period)) + 
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(limits=c("2050", "2070")) +
  theme(legend.position="none") +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Area changes for 2050 and 2070") +
  xlab("Period") +
  ylab("Future Area Change")


## Histogram for range are changes at 2070
# ggplot(data = long_data[long_data$Contraction <0 & long_data$Period == 2070,], 
#        aes(x = Contraction, fill = cut(Contraction, 10))) +
#   geom_histogram(alpha = 1, position = "identity", bins = 20) +
#   # geom_density() +
#   theme(legend.position="none") +
#   scale_fill_viridis(discrete = TRUE) +
#   ggtitle("Range contractions for 2070") +
#   ylab("Number of Species")


# note: prøv at bruge barplots: https://www.guru99.com/r-bar-chart-histogram.html

# denne figur viser sammenhæng mellem fremtidig overlap (x-akse) og hvordan range area (y-akse) ændres. 
ggplot(data = long_data, 
       aes(y = Overlap, x = Contraction, shape = Period, color = log(Distance+1))) +
  geom_point() +
  geom_line(aes(group = Species, color = log(Distance +1))) +
  scale_color_viridis(discrete = FALSE, option = "B", name = "log(Distance)") +
  theme(legend.position = c(0.85, 0.545)) +
  xlab("Area Change") + 
  ylab("Future range overlap")



# boxplot for future area changes
ggplot(data = long_data[long_data$Period != "Present",], aes(x = Elevation_Group, y = Contraction, fill = Period)) + 
  geom_boxplot(position=position_dodge(1)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("0" = "Lowland", "1"="Highland")) +
  theme(legend.position = c(0.82, 0.8)) +
  xlab("Elevation Group") +
  ylab("Area Change")


################
#  
#  Elevation changes   
#
###############test###

# how to plot subsets
# ggplot(data = long_data[long_data$Period == 2070 & long_data$Elevation > 2600,],
#        aes(x = Elevation, y = Contraction)) +
#   geom_point()

# ggplot(data = long_data, aes(x = Elevation, fill = Period)) +
#   geom_histogram(alpha = 0.6, position = "identity", bins = 40) +
#   scale_fill_viridis(discrete = TRUE) +
#   ylab("Count")


ggplot(data = long_data[long_data$Contraction >=-2,], 
       aes(x = factor(Period, level = c("Present", "2050", "2070")), 
           y = Elevation, 
           group = Species, color = Elevation_Group)) +
  geom_point(shape = 4) +
  geom_line() +
  scale_color_discrete(name = "Elevation Group", labels = c("Lowland", "Highland")) + 
  # theme(legend.position="none") +
  # scale_color_viridis(discrete = FALSE, option = "B") +
  theme(legend.position = c(0.82, 0.65)) +
  xlab("Period") +
  ylab(" Elevation (m)")


wilcox.test((long_data$Elevation[long_data$Period != "Present" & long_data$Elevation_Group == "0"]),
            (long_data$PElevation[long_data$Period != "Present" & long_data$Elevation_Group == "1"]), 
            paired = FALSE)

# data:  (long_data$Elevation[long_data$Period != "Present" & long_data$Elevation_Group == "0"])
# V = 21321, p-value < 2.2e-16
# alternative hypothesis: true location is not equal to 0




ggplot(data = long_data[long_data$Period != "Present",] , aes( x = Elevation_Group, y = Contraction)) +
  geom_boxplot() +
  # geom_line(aes(group = Species)) +
  xlab("Area change") +
  ylab("Elevation (m)")



wilcox.test((long_data$Elevation[long_data$Period != "Present" & long_data$Elevation_Group == "0"]),
            (long_data$PElevation[long_data$Period != "Present" & long_data$Elevation_Group == "1"]), 
            paired = FALSE)   



# ggplot(data = long_data[long_data$Contraction >=-2,], 
#        aes(x = factor(Period, level = c("Present", "2050", "2070")), 
#            y = Elevation, group = Species, color = Contraction)) +
# geom_boxplot()


ggplot(data = long_data, aes(x = PoleShift, y = Elevation))+
  geom_point()


#########
#  
# Range shift plots 
#  
######### 

# worldmap with range shifts  
ggplot() +
  coord_fixed(ratio = 1) +
  geom_sf(data = worldmap) +
  geom_point(data = long_data, aes(x = Centroid_X, y = Centroid_Y, 
                                   color = Elevation, shape = Period 
                                   # size = Area
  )) +
  # geom_text(data = long_data, aes(x = Elevation, y = 5, label = Elevation)) +
  geom_segment(data = long_data[long_data$Period == "2050" & long_data$Period == "Present"],  
               aes(x = long_data$Centroid_X[long_data$Period == "Present"],
                   y = long_data$Centroid_Y[long_data$Period == "Present"],
                   xend = long_data$Centroid_X[long_data$Period == "2050"],
                   yend = long_data$Centroid_Y[long_data$Period == "2050"]),
               arrow = arrow(length = unit(0.11, "cm"), type = "closed")) +
  
  geom_segment(data = long_data[long_data$Period == "2050" & long_data$Period == "2070"],  
               aes(x = long_data$Centroid_X[long_data$Period == "2050"],
                   y = long_data$Centroid_Y[long_data$Period == "2050"],
                   xend = long_data$Centroid_X[long_data$Period == "2070"],
                   yend = long_data$Centroid_Y[long_data$Period == "2070"]),
               arrow = arrow(length = unit(0.11, "cm"), type = "closed")
  ) +
  xlim(-75,-70) +
  ylim(-40,-30) +
  xlab("Longitude") +
  ylab("Latitude")  


# Scatter plot for polewards shifts and statistic test for difference between elevation groups.
ggplot(data = long_data[long_data$Period != "Present",], aes(x = PoleShift, y = Elevation, color = Period, shape = Elevation_Group)) +
  geom_point() +
  theme(legend.position = c(0.8, 0.75)) +
  geom_vline(aes(xintercept = 0)) +
  scale_shape_discrete(name = "Elevation Group", labels = c("Lowland", "Highland")) + 
  xlab("Poleshift deg.")


wilcox.test((long_data$PoleShift[long_data$Period != "Present" & long_data$Elevation_Group == "0"]),
            (long_data$Poleshift[long_data$Period != "Present" & long_data$Elevation_Group == "1"]), 
            paired = FALSE)

# data:  (long_data$PoleShift[long_data$Period != "Present" & long_data$Elevation_Group == "0"])
# V = 12833, p-value = 0.007469
# alternative hypothesis: true location is not equal to 0
# 

# plot with polewards shifts and elevation. positive y-values mean a latitudinal shif polewards, while neagtive values mean a shift equatorwards
# 
# ggplot(data = analysis_stats, aes(x = (Range_Contraction_2070), y = delta_N_total, color = Range_Contraction_2070)) +
#   geom_point() +
#   scale_color_viridis(discrete = FALSE) +
#   ylab("Total latitudinal shift") +
#   ggtitle("Polewards shift in relation with present elevation")


# Scatter plot with elevation and distance

ggplot(data = long_data[long_data$Period != "Present",], aes(x = Distance, y = Elevation, color = Period, shape = Elevation_Group)) +
  geom_point() +
  xlab("Distance (km)")

wilcox.test((long_data$Distance[long_data$Period == 2070 & long_data$Elevation_Group == "0"]),
            (long_data$Distance[long_data$Period == 2070 & long_data$Elevation_Group == "1"]), 
            paired = FALSE)

# data:  (long_data$Distance[long_data$Period != "Present" & long_data$Elevation_Group == "0"]) and (long_data$Distance[long_data$Period != "Present" & long_data$Elevation_Group == "1"])
# W = 4683, p-value = 0.001646
# alternative hypothesis: true location shift is not equal to 0



wilcox.test((long_data$Distance[long_data$Period == 2050]),
            (long_data$Distance[long_data$Period == 2070]), 
            paired = FALSE)

# data:  (long_data$Distance[long_data$Period == 2050]) and (long_data$Distance[long_data$Period == 2070])
# W = 9583, p-value = 9.41e-06


# scatterplot for poleshift and distance
ggplot(data = long_data[long_data$Period != "Present",], aes(x = PoleShift, y = log(Distance+1), color = Period)) +
  geom_point() +
  geom_line(aes(group = Species), color = "black", alpha = 0.2) +
  geom_vline(aes(xintercept = 0)) +
  theme(legend.position = c(0.8, 0.2)) +
  ylab("log(Distance) km") +
  xlab("Poleshift deg.")


# scatterplot for poleshift and distance
ggplot(data = long_data[long_data$Period != "Present",], aes(x = PoleShift_km, y = log(Distance+1), color = Period)) +
  geom_point() +
  # geom_line(aes(group = Species)) +
  geom_vline(aes(xintercept = 0)) +
  theme(legend.position = c(0.8, 0.2)) +
  ylab("log(Distance) km") +
  xlab("Poleshift km")



ggplot(data = long_data[long_data$Period != "Present",], aes(x = Period, y = log(Distance+1), fill = Period)) +
  geom_boxplot() +
  ylab("log(Distance) km")




# General scatter plot of future range shifts
# ggplot(data = long_data, aes(y = log(Distance + 1), x = Period)) +
#   geom_point(aes(alpha = 0.3), shape = 1, size = 1) +
#   geom_line(aes(group = Species, color = Shift_Accelerate)) +
#   scale_color_viridis(discrete = TRUE, option = "B") +
#   ylab("log(Distance + 1) km^2  ") +
#   xlab("Period") +
#   ggtitle("Future range shifts")

# Boxplot of range shifts  
# ggplot(data = long_data[long_data$Period == 2050,], aes(y = Distance, x = Shift_Accelerate)) +
#   geom_boxplot() +
#   # geom_point(aes(alpha = 0.3), shape = 1, size = 1) +
#   # geom_line(aes(group = Species, color = Shift_Accelerate)) +
#   # scale_color_viridis(discrete = TRUE, option = "B") +
#   ylab("log(Distance + 1) km^2  ") +
#   xlab("Period") +
#   ggtitle("Future range shifts")
# 



# plot of elevation groups in 2070 and their shifts
# creating new dataframe for this
# temp_elev_df <- long_data %>%
#   filter(Period == 2070) %>%
#   mutate(elevation_group = ifelse(Elevation > 3000, "Above 3000", "Below 3000"))
# 
# 
ggplot(data = long_data, aes(x=Elevation, y = Contraction, fill = Elevation_Group)) +
  geom_boxplot(alpha = 0.5) +
  geom_hline(yintercept = 0) +
  theme(legend.position="none") +
  scale_fill_viridis(discrete = TRUE) +
  ylab("Area change") +
  ggtitle("Area change 2070")
# 

# ggplot(data = temp_elev_df, aes(x=Elevation, y = Contraction, fill = Elevation_Group)) + 
#   geom_boxplot(alpha = 0.5) +
#   geom_hline(yintercept = 0) +
#   theme(legend.position="none") +
#   scale_fill_viridis(discrete = TRUE) +
#   ylab("Area change") +
#   ggtitle("Area change 2070")

## checking if the two elevation groups are statistically disimilar
# running wilks test since group below 3000 is not normally distributed
wilcox.test((long_data$Contraction[long_data$Elevation_Group == "0"]),
            (long_data$Contraction[long_data$Elevation_Group == "1"]), 
            paired = FALSE)

# W = 868, p-value = 0.958
# alternative hypothesis: true location shift is not equal to 0

# ggplot(data = cbind(temp_elev_df[7],temp_elev_df[9]), 
#        aes(x = Contraction, fill = elevation_group)) + 
#   geom_density(alpha = 0.5) +
#   scale_fill_viridis(discrete = TRUE) +
#   geom_vline(xintercept = 0) +
#   xlab("Area Change") +
#   ylab("Density") +
#   ggtitle("Density plot for area change in 2070.")
# 

# 
# ggplot(data = long_data, aes(x=Elevation, y = Overlap)) +
#   aes(y = Contraction, x = Overlap, shape = Elevation_Group, color = Elevation_Group) +
#   geom_point(aes(shape = Elevation_Group)) +
#   scale_shape_manual(values = c(1,4)) +
#   geom_hline(yintercept = 0, alpha = 0.5) + 
#   geom_abline(intercept = -1, slope = 1, alpha = 0.5) +
#   # geom_segment(aes(x= 0, y = -1, xend = 1, yend = 0)) +
#   # scale_color_viridis(discrete = TRUE, option = "D") +
#   ylab("Area Change") + 
#   xlab("Future range overlap") +
#   labs(shape = "Elevation Group", color = "Elevation") +
#   ggtitle("Range size changes and future overlaps")
# 
# 

##### Regular worldmap

ggplot(data = long_data[long_data$Period == 2070,]) +
  coord_fixed(ratio = 1) +
  geom_sf(data = worldmap) +
  geom_point(aes(x = Centroid_X, y = Centroid_Y, color =Elevation)) +
  # geom_text(data = long_data, aes(x = Elevation, y = 5, label = Elevation)) +
  # geom_hline(yintercept = 35) +
  scale_color_viridis(discrete = FALSE) +
  
  # xlim(72,143) +
  # ylim(25,61) +
  xlab("Longitude") +
  ylab("Latitude") 




### Extinction stuff

ggplot(data = long_data, aes(x = Elevation, y = Extinction_Overlap)) +
  geom_point() +
  geom_line(aes(group = Species))



ggplot(data = long_data[long_data$Period != "Present",], aes(x = Elevation_Group, y = Extinction_Overlap, fill = Period)) + 
  geom_boxplot(position=position_dodge(1)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("0" = "Lowland", "1"="Highland")) +
  theme(legend.position = c(0.82, 0.85)) +
  xlab("Elevation Group") +
  ylab("Extinction Risk")


wilcox.test((long_data$Extinction_Overlap[long_data$Elevation_Group == 0 & long_data$Period == 2050]),
            (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2050]), 
            paired = FALSE)

# data:  (long_data$Extinction_Overlap[long_data$Elevation_Group == 0 & long_data$Period == 2050]) and (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2050])
# W = 1227.5, p-value = 0.00811


wilcox.test((long_data$Extinction_Overlap[long_data$Elevation_Group == 0 & long_data$Period == 2070]),
            (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2070]), 
            paired = FALSE)

# data:  (long_data$Extinction_Overlap[long_data$Elevation_Group == 0 & long_data$Period == 2070]) and (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2070])
# W = 1071.5, p-value = 0.141
# 

ggplot(data = long_data[long_data$Period != "Present",], aes(x = Period, y = Extinction, fill = Period)) + 
  geom_boxplot(position=position_dodge(1)) +
  # geom_hline(yintercept = 0) +
  # scale_x_discrete(labels = c("0" = "Lowland", "1"="Highland")) +
  xlab("Period") +
  ylab("Extinction Risk")

wilcox.test((long_data$Extinction_Overlap[long_data$Elevation_Group == 0]),
            (long_data$Extinction_Overlap[long_data$Elevation_Group == 1]), 
            paired = FALSE)

# data:  (long_data$Extinction_Overlap[long_data$Elevation_Group == 0]) and (long_data$Extinction_Overlap[long_data$Elevation_Group == 1])
# W = 8877.5, p-value = 0.1377



wilcox.test((long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2050]),
            (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2070]), 
            paired = FALSE)

# data:  (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2050]) and (long_data$Extinction_Overlap[long_data$Elevation_Group == 1 & long_data$Period == 2070])
# W = 66.5, p-value = 0.007434



###########
# work place


ggplot(data = long_data[long_data$Period != "Present",], aes(x = Elevation_Group, y = log(Distance+1), fill = Period)) + 
  geom_boxplot(position=position_dodge(1)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("0" = "Lowland", "1"="Highland")) +
  xlab("Elevation Group") +
  ylab("log(Distance+1), km")




wilcox.test((long_data$Extinction[long_data$Elevation_Group == 0]),
            (long_data$Extinction[long_data$Elevation_Group == 1]), 
            paired = FALSE)

# W = 7978.5, p-value = 0.8809

wilcox.test((long_data$Extinction[long_data$Elevation_Group == 0 & long_data$Period == 2050]),
            (long_data$Extinction[long_data$Elevation_Group == 1 & long_data$Period == 2050]), 
            paired = FALSE)

# W = 1037.5, p-value = 0.2231

wilcox.test((long_data$Extinction[long_data$Elevation_Group == 0 & long_data$Period == 2070]),
            (long_data$Extinction[long_data$Elevation_Group == 1 & long_data$Period == 2070]), 
            paired = FALSE)

# W = 883, p-value = 0.9578


wilcox.test((long_data$Extinction[long_data$Period == 2050]),
            (long_data$Extinction[long_data$Period == 2070]), 
            paired = FALSE)

# W = 6407.5, p-value = 0.1398


ggplot(data = long_data[long_data$Period != "Present",], aes(x = Extinction, y= Elevation, color = Elevation_Group, shape = Elevation_Group))+
  geom_point()



ggplot(data = long_data, aes(x = Overlap, y = log(Distance+1), color = Elevation_Group)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE)


ggplot(data = long_data[long_data$Period != "Present" & long_data$Contraction < 0,], aes(x = Contraction, color = Period)) +
  
  
  # x axis treated as continuous variable
  # df2$dose <- as.numeric(as.vector(df2$dose))
  # ggplot(data=df2, aes(x=dose, y=len, fill=supp)) +
  #   geom_bar(stat="identity", position=position_dodge())+
  #   scale_fill_brewer(palette="Paired")+
  #   theme_minimal()
  
  
  

ggplot(data=long_data[long_data$Period != "Present" & long_data$Contraction < 0,], aes(x = Period, y = Contraction)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")


# ggplot(data = long_data[long_data$Period != "Present",], aes(x = Contraction, color = Period)) +
# geom_bar(aes(fill = Period, stat = Period), position=position_dodge()) +
#   scale_x_continuous(limits = c(-1,0))


#
ggplot(data = long_data[long_data$Period != "Present",], aes(x = log(Distance+1), y = PoleShift, color = Period)) +
  geom_point(aes(fill = Period))


# ggplot(data = long_data[long_data$Contraction >=-2,], aes(x = factor(Period, level = c("Present", "2050", "2070")), 
#                                                           y = Elevation, group = Species, color = log(Distance+1))) +
#   geom_point() +
#   geom_line() +
#   # theme(legend.position="none") +
#   scale_color_viridis(discrete = FALSE, option = "B") +
#   xlab("Time Periods")
# 

ggplot(data = long_data[long_data$Period != "Present",], aes(x = Centroid_Y, y = Contraction)) + 
  geom_point(aes( color = Elevation_Group)) +
  # geom_line(aes(group = Species)) +
  # geom_boxplot(position=position_dodge(1)) +
  # geom_hline(yintercept = 0) +
  # scale_x_discrete(labels = c("0" = "Lowland", "1"="Highland")) +
  xlab("Latitude") +
  ylab("Area Change")



ggplot(data = long_data[long_data$Period != "Present",], 
       aes(x = log(Area+1), fill = Period)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 40) +
  scale_fill_viridis(discrete = TRUE) +
  xlab("log-Area") + 
  ylab("Count")




###  

## Plots for elevation in relation to distance

plot(analysis_stats$Present_Mean_Elevation, analysis_stats$dists_2050, 
     xlab = "Mean Elevation Present",
     ylab = "Dists 2050")

plot(analysis_stats$Present_Area, analysis_stats$dists_2050, 
     xlab = "Present Range Area(Km^2)",
     ylab = "Dists 2050 (Km)")  


plot(analysis_stats$Mean_Elevation_2070, analysis_stats$dists_2070, 
     xlab = "Mean Elevation 2070",
     ylab = "Dists 2070")



# making histogram for elevation developments

ggplot(data = analysis_stats) +
  geom_histogram(aes(x = Present_Mean_Elevation), fill = "red", binwidth = 50) +
  geom_histogram(aes(x = Mean_Elevation_2050), alpha = 0.5, binwidth = 50) +
  geom_histogram(aes(x = Mean_Elevation_2070), alpha = 0.5, binwidth = 50) +
  xlab("Elevations")  



# scatter plot with present day elevation and future range shifts

ggplot(data = analysis_stats, aes(x = Present_Mean_Elevation, y = dists_2050)) +
  geom_point() +
  geom_smooth(method = lm)


# scatter plot with present day area and future range shifts
ggplot(data = analysis_stats, aes(x = Present_Area, y = dists_2050)) +
  geom_point() +
  geom_smooth(method = lm)


## Making scatterplot with sum of distances and elevation in 2070
ggplot(data = analysis_stats, aes(x = (dists_2050 + dists_2070), 
                                  y = Mean_Elevation_2070)) +
  geom_point() +
  geom_smooth(method = lm)







# ggplot(data = analysis_stats) + 
#   coord_fixed(ratio = 1) + 
#   wm +
#   geom_point(aes(x = `Centroid_x Present`, y = `Centroid_y Present`),
#              size = 1) +
#   geom_point(aes(x = `Centroid_x 2050`, y = `Centroid_y 2050`),
#              size = 1, color = "red") +
#   geom_point(aes(x = `Centroid_x 2070`, y = `Centroid_y 2070`),
#              size = 1, color = "yellow") +
#   xlim(-20,  20) +
#   ylim(30,  60) +
#   theme_bw()  


# plot with centroid movements and point size is area size  
ggplot() +
  coord_fixed(ratio = 1) +
  wm +
  geom_point(data = analysis_stats, 
             aes(x = Centroid_x_Present,
                 y = Centroid_y_Present,
                 size = Present_Area)) +
  
  geom_point(data = analysis_stats,
             aes(x = Centroid_x_2050,
                 y = Centroid_y_2050)) +
  
  geom_point(data = analysis_stats,
             aes(x = Centroid_x_2070,
                 y = Centroid_y_2070,
                 size = Area_2070)) +
  # scale_color_viridis() +
  
  geom_segment(data = analysis_stats,
               aes(x = Centroid_x_Present,
                   y = Centroid_y_Present,
                   xend = Centroid_x_2050,
                   yend = Centroid_y_2050,
                   colour = "red"),
               arrow = arrow(angle = 30,
                             length = unit(0.07, "inches"),
                             ends = "last",
                             type = "closed")) +
  
  geom_segment(data = analysis_stats,
               aes(x = Centroid_x_2050,
                   y = Centroid_y_2050,
                   xend = Centroid_x_2070,
                   yend = Centroid_y_2070,
                   colour = "black"),
               arrow = arrow()) +
  
  # xlim(-15,20)+
  # ylim(30,65)+
  labs(title = "not log transformed") + 
  theme(legend.position = "none")


