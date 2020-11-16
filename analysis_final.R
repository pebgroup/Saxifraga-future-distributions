# Script to analyze y model outputs

## Libraries

library(raster)
library(tidyverse)


# import of elevation raster

elevation <- raster("wc2.1_2.5m_elev//wc2.1_2.5m_elev.tif", xmn=-180, xmx=180, ymn=-60, ymx=90)

## Part where I import files from my model directory


folders <- list.files("data//models")
length(folders)

# creating empty vectors
ranges_present <- vector()
ranges_2050 <- vector()
ranges_2070 <- vector()

centroid_present_x <- c()
centroid_present_y <- c()
centroid_2050_x <- c()
centroid_2050_y <- c()
centroid_2070_x <- c()
centroid_2070_y <- c()

elevation_present <- vector()
elevation_2050 <- vector()
elevation_2070 <- vector()
elevation_present_sd <- vector()
elevation_2050_sd <- vector()
elevation_2070_sd <- vector()


overlaps_2050 <- vector()
overlaps_2070 <- vector()

change_2050 <- vector()
change_2070 <- vector()

AUCs <- vector()
AUCs_sd <- vector()




for(i in 1:length(folders)){
  files <- list.files(paste0("data//models//", folders[i]))
  species_name <- folders[i]
  files_list <- list()
  # In the event that some species folders are empty, I include this small check to skip empty folders.    
  if(length(files)==0){
    print(paste(i,folders[i], "is empty. Skipping."))
  } else {
    print(paste(i,"Analyzing", species_name))
    
    
    for(j in 1:length(files)){
      files_list[[j]] <- readRDS(paste0("data//models//", folders[i],"//", files[j]))
      species_name <- folders[i]
    }
    names(files_list) <- files
    
    suitability_present <- files_list$suitability_present.rds > 0.7
    suitability_2050 <- files_list$suitability_2050.rds > 0.7
    suitability_2070 <- files_list$suitability_2070.rds > 0.7
    stats <- files_list$stats.rds
    
    
    ## area analysis
    # calculate areas of present and future suitabilities
    # note: the areas are not entirely accurate - i think this is due to which map projection I have chosen.
    # i will recalculate to dimension-less values later on.
    
    area_present <- suitability_present
    area_present[area_present <1] <- NA
    area_present <- area(area_present, na.rm = TRUE, weights = FALSE)
    area_present <- area_present[!is.na(area_present)]
    area_present <-length(area_present)*median(area_present) #output in square km
    
    area_2050 <- suitability_2050
    area_2050[area_2050 <1 ] <- NA
    area_2050 <- area(area_2050, na.rm = TRUE, weights = FALSE)
    area_2050 <- area_2050[!is.na(area_2050)]
    area_2050 <- length(area_2050)*median(area_2050) #output in square km
    
    area_2070 <- suitability_2070
    area_2070[area_2070 <1 ] <- NA
    area_2070 <- area(area_2070, na.rm = TRUE, weights = FALSE)
    area_2070 <- area_2070[!is.na(area_2070)]
    area_2070 <- length(area_2070)*median(area_2070) #output in square km
    
    
    # calculate center of the suitabilities
    
    center_present <- colMeans(xyFromCell(suitability_present, which(suitability_present[]==1)))
    center_2050 <- colMeans(xyFromCell(suitability_2050, which(suitability_2050[]==1)))
    center_2070 <- colMeans(xyFromCell(suitability_2070, which(suitability_2070[]==1)))
    
    # calculate elevation ranges and mean
    
    elevation_extent_present <- extent(suitability_present)
    elevation_present_crop <- raster::crop(elevation, elevation_extent_present)
    elevation_range_present <- elevation_present_crop[suitability_present == 1,]
    
    # Observations falling into the ocean will be given below sealevel elevations. I rewrite those values to 0.
    elevation_range_present[elevation_range_present <0] <- 0
    
    sd(elevation_range_present, na.rm = TRUE)
    mean(elevation_range_present, na.rm = TRUE)
    
    
    elevation_extent_2050 <- extent(suitability_2050)
    elevation_2050_crop <- raster::crop(elevation, elevation_extent_2050)
    elevation_range_2050 <- elevation_2050_crop[suitability_2050 == 1,]
    elevation_range_2050[elevation_range_2050 <0] <- 0
    
    sd(elevation_range_2050, na.rm = TRUE)
    mean(elevation_range_2050, na.rm = TRUE)
    
    
    
    elevation_extent_2070 <- extent(suitability_2070)
    elevation_2070_crop <- raster::crop(elevation, elevation_extent_2070)
    elevation_range_2070 <- elevation_2070_crop[suitability_2070 == 1,]
    elevation_range_2070[elevation_range_2070 <0] <- 0
    
    # sd(elevation_range_2070, na.rm = TRUE)
    # mean(elevation_range_2070, na.rm = TRUE)
    
    
    
    ## determine overlap
    # note: rasters require the same extent to determine overlap.
    overlap_2050 <- overlay(suitability_2050, suitability_present, fun = sum)
    overlap_2050 <- overlap_2050 == 2
    
    overlap_2070 <- overlay(suitability_2070, suitability_present, fun = sum)
    overlap_2070 <- overlap_2070 == 2
    
    
    # calculate area of overlaps
    
    overlap_area_2050 <- overlap_2050
    overlap_area_2050[overlap_area_2050 <1 ] <- NA
    overlap_area_2050 <- area(overlap_area_2050, na.rm = TRUE, weights = FALSE)
    overlap_area_2050 <- overlap_area_2050[!is.na(overlap_area_2050)]
    overlap_area_2050 <- length(overlap_area_2050)*median(overlap_area_2050) #output in square km
    if(is.empty(overlap_area_2050) == TRUE){overlap_area_2050 <- 0}
    
    overlap_area_2070 <- overlap_2070
    overlap_area_2070[overlap_area_2070 <1 ] <- NA
    overlap_area_2070 <- area(overlap_area_2070, na.rm = TRUE, weights = FALSE)
    overlap_area_2070 <- overlap_area_2070[!is.na(overlap_area_2070)]
    overlap_area_2070 <- length(overlap_area_2070)*median(overlap_area_2070) #output in square km
    if(is.empty(overlap_area_2070) == TRUE){overlap_area_2070 <- 0}
    
    
    # Calculating range spatial changes
    
    if(is.empty(area_2050) == TRUE){change_2050[i] <- -1} else {
      change_2050[i] <- round(((area_2050 - area_present) / area_present), digits = 2)} 
    
    if(is.empty(area_2070) == TRUE){change_2070[i] <- -1} else {
      change_2070[i] <- round(((area_2070 - area_present) / area_present), digits = 2)} 
    
    
    
    
    
    
    # appending data
    
    ranges_present[i] <- round(area_present, digits = 2)
    
    if(is.empty(area_2050) == TRUE){ranges_2050[i] <- 0} else {
      ranges_2050[i] <- round(area_2050, digits = 2)}
    
    if(is.empty(area_2070) == TRUE){ranges_2070[i] <- 0} else {
      ranges_2070[i] <- round(area_2070, digits = 2)}
    
    centroid_present_x[i] <- center_present[1]
    centroid_present_y[i] <- center_present[2]
    
    centroid_2050_x[i] <- center_2050[1]
    centroid_2050_y[i] <- center_2050[2]
    
    centroid_2070_x[i] <- center_2070[1]
    centroid_2070_y[i] <- center_2070[2]
    
    
    elevation_present[i] <- round(mean(elevation_range_present, na.rm = TRUE), digits = 2)
    elevation_2050[i] <- round(mean(elevation_range_2050, na.rm = TRUE), digits = 2)
    elevation_2070[i] <- round(mean(elevation_range_2070, na.rm = TRUE), digits = 2)
    
    
    elevation_present_sd[i] <- round(sd(elevation_range_present, na.rm = TRUE), digits = 2)
    elevation_2050_sd[i] <- round(sd(elevation_range_2050, na.rm = TRUE), digits = 2)
    elevation_2070_sd[i] <- round(sd(elevation_range_2070, na.rm = TRUE), digits = 2)
    
    overlaps_2050[i] <- overlap_area_2050
    overlaps_2070[i] <- overlap_area_2070
    
    AUCs[i] <- round(mean(stats$AUCs, na.rm = TRUE), digits = 2)
    AUCs_sd[i] <- round(sd(stats$AUCs, na.rm = TRUE), digits = 2)
  }
}


# Combining the results into a dataframe 
analysis_stats <- tibble(folders,                           
                         ranges_present, ranges_2050, ranges_2070, 
                         centroid_present_x, centroid_present_y,
                         centroid_2050_x, centroid_2050_y,
                         centroid_2070_x, centroid_2070_y,
                         elevation_present, elevation_present_sd,
                         elevation_2050, elevation_2050_sd,
                         elevation_2070, elevation_2070_sd,
                         overlaps_2050, overlaps_2070,
                         change_2050, change_2070,
                         AUCs, AUCs_sd
)


names(analysis_stats) <- c("Species", 
                           "Present_Area", "Area_2050", "Area_2070", 
                           "Centroid_x_Present", "Centroid_y_Present", 
                           "Centroid_x_2050", "Centroid_y_2050", 
                           "Centroid_x_2070", "Centroid_y_2070",
                           "Present_Mean_Elevation", "Present_Mean_Elevation_Standard_Deviation",
                           "Mean_Elevation_2050", "Mean_Elevation_2050_Standard_Deviation",
                           "Mean_Elevation_2070", "Mean_Elevation_2070_Standard Deviation",
                           "Area_overlap_2050", "Area_overlap_2070", 
                           "Range_Contraction_2050", "Range_Contraction_2070",
                           "AUC_mean", "AUC_Standard_Deviation"
)

####  
saveRDS(analysis_stats, file = "data//analysis//analysis_stats.rds")
write.csv(analysis_stats, file = "data//analysis//analysis_stats.csv")


#############  

plot(analysis_stats$`Centroid_x Present`,analysis_stats$`Centroid_y Present`)
points(analysis_stats$`Centroid_x 2050`, analysis_stats$`Centroid_y 2050`, col = "red")
points(analysis_stats$`Centroid_x 2070`, analysis_stats$`Centroid_y 2070`, col = "yellow")

#####  
stats_fragilis <- readRDS(paste0("data//models//Saxifraga fragilis//stats.rds"))
stats_granulata <- readRDS(paste0("data//models//Saxifraga granulata//stats.rds"))
stats_hypnoides <- readRDS(paste0("data//models//Saxifraga hypnoides//stats.rds"))
stats_strigosa <- readRDS(paste0("data//models//Saxifraga strigosa//stats.rds"))

sd(stats_fragilis$AUCs)
sd(stats_granulata$AUCs)
sd(stats_hypnoides$AUCs)
sd(stats_strigosa$AUCs)

test <- files_list$stats.rds


round(sd(readRDS(paste0("data//models//Saxifraga strigosa//stats.rds"))$AUCs), digits = 4)

# files_list$stats.rds$AUCs

#
#List file subdirectories
folders<- list.files(path = "yourfolder")
#
#Get all files...
files <- rep(NA,0)
for(i in c(1:length(folders))){
  files.i <- list.files(path = noquote(paste("yourfolder/",folders[i], "/", sep = ""))) 
  n <- length(files.i)
  files.i <- paste(folders[i], files.i, sep = "/")
  files <- c(files, files.i)
} 
# 
#
#Read first data file (& add file name as separate column)
T1 <- read.delim(paste("yourfolder/", files[1], sep = ""), sep = "", header=TRUE)
T1 <- cbind(T1, "FileName" = files[1])