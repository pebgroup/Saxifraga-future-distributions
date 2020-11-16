## Species distribution modeling script and future model

# Libraries
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(dismo)
library(spatialEco)
library(biogeo)

library(blockCV)
library(maxnet)

library(precrec)
library(ecospat)
library(usdm)

library(retry)




# Inloading distribution matrix, and species list, which are necessary for cleaning coordinates based on known species ranges.
comm <- readRDS("data//Shapefiles//comm.rds")

# Inloading shapefile for the world. This is needed for range filtering.
world <- st_read("data//Shapefiles//level3.shp")  

# Inloading climate data 
climate <- getData(name = "worldclim", 
                   var = "bio", 
                   res = 2.5)

# Inloading future predicted climate scenarios here.
# I choose two scenarios based on the following:
# Model: The Norwegian Earth System Model, NorESM1-M
# Emission scenario: 4.5
# Years are for 2050 and 2070  

future_climate_2050 <- getData('CMIP5', var='bio', res=2.5, rcp=45, model='NO', year=50)
names(future_climate_2050) <- names(climate)

future_climate_2070 <- getData('CMIP5', var='bio', res=2.5, rcp=45, model='NO', year=70)
names(future_climate_2070) <- names(climate)




# Picking the variables I want to use for my modeling. 
# These are: bio5; max temperature of warmest month, bio6; min temperature of coldest month, bio13; precipitation of wettest month, bio14; precipitation of driest month
climate <- climate[[c("bio1", "bio5", "bio7", "bio12", "bio18")]]
future_climate_2050 <- future_climate_2050[[c("bio1", "bio5", "bio7", "bio12", "bio18")]]
future_climate_2070 <- future_climate_2070[[c("bio1", "bio5", "bio7", "bio12", "bio18")]]

# Running variance inflation on my chosen variables to see if there is collinerarity. Variables scoring over 5 have collinerarity.  
vif(climate)

# Variables        VIF
# 1      bio1 262.710074
# 2      bio5 140.908302
# 3      bio7  60.517714
# 4     bio12   4.456075
# 5     bio18   3.696752

# Looks ´like there are issues. I remove bio1, since it scores highest, and run the test again

climate <- climate[[c("bio5", "bio7", "bio12", "bio18")]]

vif(climate)
# Variables      VIF
# 1      bio5 1.175290
# 2      bio7 1.914092
# 3     bio12 3.956519
# 4     bio18 2.900679


# Okay, now the remaining variables all score below 5, so I proceed.

# First I rewrite the future climate variables to reflect the changes I ust did.

future_climate_2050 <- future_climate_2050[[c("bio5", "bio7", "bio12", "bio18")]]
future_climate_2070 <- future_climate_2070[[c("bio5", "bio7", "bio12", "bio18")]]


# And just to make double sure the future climate variables do not share collinerarity, I run a vif test on them too.
# vif(future_climate)
# Variables      VIF
# 1      bio5 1.107382
# 2      bio7 1.718191
# 3     bio12 3.430377
# 4     bio18 2.561970

# Okay everything checks out.  

# Making sure all climate models have the same crs

crs(future_climate_2050) <- crs(climate)
crs(future_climate_2070) <- crs(climate)

# First part of the script
# Here data is inloaded and prepared for the modeling.
# The data has been previously cleaned, and all presence points are assumed to have been associated with environmental background data.

dat <- read.csv("data//Species//all_saxifrage_occurrences_cleaned_12_10.csv", header = TRUE)

# work, go with dataset of 5 species
# granulata: 58863
# hypnoides: 3824
# fragilis: 491
# strigosa: 50
# glabella: 10

# dat <- dat[dat$Species == "Saxifraga capitata",]
#            dat$Species == "Saxifraga hypnoides" |
#            dat$Species == "Saxifraga fragilis"  |
#            dat$Species == "Saxifraga strigosa"  |
#            dat$Species == "Saxifraga glabella",]
# 

### testing out species that had corrupted suitability_2050 during previous modeling
# 
# dat <- dat[dat$Species == "Saxifraga aizoides" | 
#            dat$Species == "Saxifraga cernua" |
#            dat$Species == "Saxifraga hirculus"  |
#            dat$Species == "Saxifraga oppositifolia"  |
#            dat$Species == "Saxifraga rivularis",]


# Thinning data using the duplicatesexclude function from the biogeo library
# I also use the function missingvalsexclude to remove any observations that do not fall inside the climate raster.
# These should all have been removed in the cleaning script - I do it again to make sure.  
dat <- dat %>%
  transmute(ID = seq(nrow(dat)),
            Species, 
            decimallongitude, 
            decimallatitude,
            coordinateUncertaintyInMeters,
            elevationData,
            year,
            basisOfRecord, 
            countryCode,
            taxonRank,
            plantNameID,
            gbifID)


dat <- keepmainfields(dat,
                      ID = "ID",
                      Species = "Species",
                      x = "decimallongitude",
                      y = "decimallatitude",
                      others = c("elevationData", "year", "plantNameID"))

dat <- duplicatesexclude(dat, res = 2.5)

dat <- missingvalsexclude(climate[[1]], dat)

dat <- dat %>%
  filter(Exclude == 0)

# If any species drops below n=10, I remove them from the dataset.
species_count <- speciescount(dat, orderby = "Species")

accepted_species <- species_count %>%
  filter(ntot > 10)

####
# costum filter bit here for diagnostic purposes or divying up the workload  


# only modelling the first 20 species  
# accepted_species <- accepted_species[-(1:78),]


##### end of costum filter bit  

dat <- dat %>%
  filter(Species %in% accepted_species$Species) %>%
  transmute(species = Species,
            decimallongitude = x,
            decimallatitude = y,
            year,
            elevationData,
            plantNameID) 




#############################
#  
# Recreate the dataframe as a list, with each list element being a seperate species.
#  
################################  
data_list <- list()
species_names <- unique(dat$species)

# sorting speciesnames alpabetically
species_names <- sort(species_names)

for(i in 1:length(species_names)){
  data_list[[i]] <- subset(dat, dat$species == species_names[i])
  names(data_list)[i] <- species_names[i]
}

# Create a duplicate list where the elements are turned into sf format.  

data_list_sf <- data_list

for (i in 1:length(species_names)){
  data_list_sf[[i]] <- data_list_sf[[i]] %>%
    st_as_sf(coords = c("decimallongitude", "decimallatitude")
    ) %>%
    as_Spatial()
  crs(data_list_sf[[i]]) <- crs(climate)
  
}


###########  
#  
# Full model script here.
#  
##################  


# Beginning by creating absence points for each species.  

plant_nameID <- vector()
# distribution_list <- list()
# climate_mask_list <- list()
# climate_data_list <- list()
# models_list <- list()


for(i in 1:length(species_names)){
  # I want to save the models for each species in their own directory.
  # So first I make an if loop that checks if a directory does not exist. If it does not exist the loop then creates the directory.  
  if(!dir.exists(paste0("data//models//",species_names[i]))){
    dir.create(paste0("data//models//",species_names[i]))
    print(paste0("Directory for ", species_names[i], " created."))}
  
  start_time <- Sys.time() # Adding timing function to see how long it takes for each species.
  print(paste(i,"Beginning modeling of", species_names[i]))
  
  #  Pulling plant name ID from dataset in preparation for looking up species range shapefiles.       
  plant_nameID <- unique(data_list_sf[[i]]$plantNameID)
  country_list <- comm[, plant_nameID]
  country_names <- names(country_list[country_list == 1])
  
  #  Subsetting the world shapefile to get the ranges for the species.    
  world_subset <- world %>%
    dplyr::filter(LEVEL_3_CO %in% country_names)
  
  #  Creating future ranges by buffering the current range.
  #  The buffer uses degrees as units. 200km is around 1.8 degrees at the equator.   
  future_buffer <- world_subset[1] %>%
    st_buffer(dist = 1.8)      
  
  future_buffer_extent <- raster::extent(x = future_buffer)
  
  print(paste("Creating climate mask for", species_names[i]))
  
  # Masking out the climate raster in accordance with the species range shapefile. 
  # Here I choose to assign the extent from the future scenario my present range. This way the result plots wlll have the same dimensions.    
  world_subset_extent  <- future_buffer_extent
  climate_crop <- raster::crop(x = climate, y = world_subset_extent)
  climate_mask <- raster::mask(x = climate_crop, mask = world_subset) 
  
  print(paste0("Done creating climate mask for ", species_names[i]))
  
  
  # Creating raster mask for future climate scenarios. 
  future_climate_2050_crop <- raster::crop(future_climate_2050, future_buffer_extent)
  future_climate_2070_crop <- raster::crop(future_climate_2070, future_buffer_extent)
  
  future_climate_2050_mask <- raster::mask(x = future_climate_2050_crop, mask = future_buffer)
  future_climate_2070_mask <- raster::mask(x = future_climate_2070_crop, mask = future_buffer)
  
  # climate_mask_list[[i]] <- climate_mask
  # names(climate_mask_list)[i] <- species_names[i]
  
  ## Creating absence points  
  print(paste(i,"Assigning background points for", species_names[i]))
  absence_points <- randomPoints(mask = climate_mask[[1]], 
                                 n = nrow(data_list_sf[[i]]), 
                                 p = data_list_sf[[i]], 
                                 lonlatCorrection = TRUE,
                                 warn = 0) %>%
    as.data.frame() %>% 
    dplyr::select(x = x, 
                  y = y) %>%
    add_column(species = 0)
  
  
  
  
  ## Creating dataframe with presence points and species set to 1
  presence_points <- data_list[[i]] %>%
    transmute(x = decimallongitude,
              y = decimallatitude,
              species = 1
    )
  
  ## Merging absencepoints with presence points into same dataframe.                         
  all_points_data <-  bind_rows(presence_points, absence_points) %>%
    transmute(ID = seq.int((nrow(presence_points) + nrow(absence_points))),
              decimallongitude = x,
              decimallatitude = y,
              species)     
  
  ## Occasionally one or more points will fall outside the climate raster, and therefore introduce NA values in the climate data later on.
  ## I run biogeo's function to get rid of those points.       
  all_points_data <- keepmainfields(all_points_data,
                                    ID = "ID",
                                    Species = "species",
                                    x = "decimallongitude",
                                    y = "decimallatitude")
  
  all_points_data <- missingvalsexclude(climate[[1]], all_points_data)
  all_points_data <- all_points_data %>%
    filter(Exclude == 0) %>%
    transmute(species = Species,
              decimallongitude = x,
              decimallatitude = y,
              ID) %>%
    st_as_sf(x = . , coords = c("decimallongitude", "decimallatitude")) %>%
    as_Spatial()
  
  crs(all_points_data) <- crs(climate)    
  
  # distribution_list[[i]] <- all_points_data
  # names(distribution_list)[i] <- species_names[i]
  
  
  
  ### Creating climate data
  
  climate_data <- raster::extract(climate, all_points_data, method = 'simple', df = TRUE)
  climate_data <- climate_data[-1]
  
  # climate_data_list[[i]] <- climate_data
  # names(climate_data_list)[i] <- species_names[i]
  
  
  
  ######## 
  #
  #  Creating full model.
  #
  ########
  
  full_model <- maxnet(
    p = all_points_data$species,
    data = climate_data,
    f = maxnet.formula(p = all_points_data$species, 
                       data = climate_data, 
                       classes = 'lqp'))
  
  
  # Creating suitability raster  
  final_suitability <- predict(object = climate_crop, 
                               model = full_model,
                               type = "cloglog")
  
  
  # Creating future predicted suitability plots for both scenarios.
  future_suitability_2050 <- predict(object = future_climate_2050_crop,
                                     model = full_model,
                                     type = "cloglog")    
  
  
  future_suitability_2070 <- predict(object = future_climate_2070_crop,
                                     model = full_model,
                                     type = "cloglog")   
  
  
  ###
  # Performing model evaluation
  ### 
  
  # Occasionally spatial_blocks will throw an error: 
  #     data.frame(blocks = seq_len(nrowBlocks), folds = systematicNum(subBlocks,  : arguments imply differing number of rows: 149, 150   
  # I don't understand in detail why or when it happens, but it seems to go away if I retry enough times. 
  # Therefore I am using the retry function from the retry library to re-run spatial_blocks untill it works. 
  # I have set the maximum number of retries to 100. Hopefully this should suffice.  
  
  spatial_blocks <- retry(
    blockCV::spatialBlock(speciesData = all_points_data, # presence-background data
                          species = "species",
                          rasterLayer = climate_mask,
                          rows = 19,
                          cols = 20,
                          k = 4,
                          selection = "random",
                          iteration = 10,
                          numLimit = 2,
                          biomod2Format = FALSE),
    when = "replacement",
    upon = "error",
    until = ~TRUE,
    max_tries = 100)
  
  AUCs <- vector()
  MSSs <- vector()
  folds <- spatial_blocks$folds
  
  # species_temp <- all_points_data$species
  # climate_temp <- climate_data
  
  for(k in 1:length(folds)){
    trainSet <- unlist(folds[[k]][1]) # extract the training set indices
    testSet <- unlist(folds[[k]][2]) # extract the testing set indices
    
    ## Creating model with maxent
    points_train <- all_points_data$species[trainSet]
    climate_train <- climate_data[trainSet, ]
    
    if(sum(points_train) < 2){next} else {
      
      mx <- maxnet(points_train, 
                   climate_train, 
                   maxnet.formula(points_train, 
                                  climate_train, 
                                  classes = "lqp"))
      
      
      testTable <- all_points_data[testSet, ] # a table for testing predictions and reference data
      testTable$pred <- predict(mx, 
                                climate_data[testSet, ], 
                                type = "cloglog") # predict the test set
      
      # evaluation statistics preparation
      
      ##  Calculate MSS
      
      prediction <- predict(mx,
                            climate_data,
                            type = 'cloglog')
      
      
      eval <- dismo::evaluate(p = prediction[all_points_data$species == 1],
                              a = prediction[all_points_data$species == 0])
      MSS <- eval@t[which.max(eval@TPR + eval@TNR)]
      
      
      # calculate area under the curve (AUC) of the receiver operating curve using the precrec package
      
      if(sum(testTable$species) == 0){
        auc <- NA} else{
          precrec_obj <- precrec::evalmod(scores = testTable$pred,
                                          labels = testTable$species)
          auc <- auc(precrec_obj)[1, 4]
        }
      
      # # from the iterations
      AUCs[k] <- as.numeric(auc)
      MSSs[k] <- as.numeric(MSS)
      
    }}
  stats <- as.data.frame(cbind(AUCs, MSS))
  
  # stats_list[[i]] <- stats
  # names(stats_list)[i] <- species_names[i]
  
  
  
  #### Done with model evaluation  
  
  
  
  
  
  
  
  # Saving data.
  
  saveRDS(all_points_data, file = paste0("data//models//",species_names[i], "//", "distribution.rds"))
  saveRDS(climate_data, file = paste0("data//models//",species_names[i], "//", "climate_data.rds"))
  saveRDS(climate_crop, file = paste0("data//models//",species_names[i], "//", "climate_crop.rds"))
  
  saveRDS(future_climate_2050, file = paste0("data//models//",species_names[i], "//", "climate_2050.rds"))
  saveRDS(future_climate_2070, file = paste0("data//models//",species_names[i], "//", "climate_2070.rds"))
  saveRDS(future_climate_2050_crop, file = paste0("data//models//",species_names[i], "//", "climate_crop_2050.rds"))
  saveRDS(future_climate_2070_crop, file = paste0("data//models//",species_names[i], "//", "climate_crop_2070.rds"))
  
  saveRDS(full_model, file = paste0("data//models//",species_names[i], "//", "SDM.rds"))
  
  saveRDS(final_suitability, file = paste0("data//models//",species_names[i], "//", "suitability_present.rds"))
  saveRDS(future_suitability_2050, file = paste0("data//models//",species_names[i], "//", "suitability_2050.rds"))
  saveRDS(future_suitability_2070, file = paste0("data//models//",species_names[i], "//", "suitability_2070.rds"))
  
  saveRDS(stats, file = paste0("data//models//",species_names[i], "//", "stats.rds"))
  saveRDS(world_subset, file = paste0("data//models//",species_names[i], "//", "range_shapefile.rds"))
  
  # Removing data before next iteration.   
  rm(absence_points, 
     all_points_data,
     climate_crop,
     climate_data,
     climate_mask,
     eval,
     final_suitability,
     folds,
     full_model,
     future_buffer,
     future_buffer_extent,
     future_climate_2050_crop,
     future_climate_2050_mask,
     future_climate_2070_crop,
     future_climate_2070_mask,
     future_suitability_2050,
     future_suitability_2070,
     mx,
     precrec_obj,
     prediction,
     presence_points,
     stats,
     testTable,
     world_subset,
     world_subset_extent)
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time   
  print(paste0("Done modeling for ", species_names[i]))
  print(time_taken)
  
} 

# end of sdm script
