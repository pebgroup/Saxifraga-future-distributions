# cleaning script, final edition.

# libraries
library(tidyverse)
library(Taxonstand)
library(raster)
library(countrycode)
library(sp)
library(rworldmap)  
library(biogeo)
library(sf)


# inloading data

species <- readRDS("data//Shapefiles//WCSP_clean.apg.rds")
# Filtering for Saxifraga only.  
sax_species <- species[species$genus == "Saxifraga",] 


dat <- read.csv("data//gbif_raw_data//0059292-200613084148143.csv", sep = "\t", header = TRUE, quote = "") #gbif filter: taxonkey: Saxifraga L.



# Begin by filtering out unwanted columns, and observations with no coordinates
dat <- dat %>% 
  transmute(
    ID = seq.int(nrow(dat)),
    species, 
    countryCode,
    decimallongitude = as.numeric(decimalLongitude), 
    decimallatitude = as.numeric(decimalLatitude),
    coordinateUncertaintyInMeters,
    year,
    gbifID, 
    elevation,
    elevationAccuracy, 
    basisOfRecord,
    taxonRank,
    speciesKey,
    taxonKey,
    genus
  )  %>%
  filter(
    !is.na(decimallongitude),
    !is.na(decimallatitude),
    basisOfRecord == "MACHINE_OBSERVATION" |
      basisOfRecord == "HUMAN_OBSERVATION" |
      basisOfRecord == "OBSERVATION" |
      basisOfRecord == "PRESERVED_SPECIMEN" |
      basisOfRecord == "LIVING_SPECIMEN",
    taxonRank == "SPECIES" |
      taxonRank == "SUBSPECIES" |
      taxonRank == "VARIETY" |
      taxonRank == "FORM"
  )


# With the filtering done, it is time to do a few diagnostics on the data, 
# such as investigating how many species are in the data, and if their names are correct.  


# There is a small naming issue with 1 species that has to be handled manually, otherwise the following functions will not perform accurately.
# dat$species[dat$species == "Saxifraga Ã— arendsii Ã— granulata"] <- "Saxifraga × arendsii × granulata"  

### Checking that all species have been properly taxonomically resolved, by running species names through The Plant List.

# Creating species list  
# species_list <- unique(dat$species)



# 
#   
#     
# # # Looking up names on The Plant List  
# #   species_list_resolved <- TPL(species_list)
# # 
# # # Making new species list with the verified species names.  
# #   verified_species <- species_list_resolved$New.Species
# #   
# # Making a loop to look up the verified species in the sax_species dataframe and return the correct plant name IDs.  
# # Defining empty containers for the loop
#   plant_name_ID_temp_frame <- data.frame()
#   plant_name_ID_list <- vector()
#   
# #Running the loop.
#   for (i in 1:(length(verified_species))){
#     plant_name_ID_temp_frame <- sax_species %>%
#       filter(species == verified_species[i] &
#                taxon_rank == "Species" &
#                taxon_status == "Accepted")
#     
#     plant_name_ID_temp_frame <- plant_name_ID_temp_frame$plant_name_id
#     
#     plant_name_ID_temp_frame <- as.character(plant_name_ID_temp_frame) 
#     if (is_empty(plant_name_ID_temp_frame) == TRUE ){
#       plant_name_ID_list[i] <- NA
#       print(paste(plant_name_ID_list[i], "rejected."))}  else {
#       plant_name_ID_list[i] <- plant_name_ID_temp_frame
#       print(paste(plant_name_ID_list[i], "Accepted."))
#       
#     }}
#   
# 
# # Creating a new dataframe with species names, plant ids, and hybrid information.
#   verified_species_frame <- as.data.frame(cbind(verified_species, plant_name_ID_list, species_list_resolved$Species, species_list_resolved$New.Hybrid.marker, species_list_resolved$New.Genus))
#   names(verified_species_frame) <- c("verifiedSpeciesName", "plantNameID", "speciesName", "hybrid", "verifiedGenus")
#   
# # In case any species were duplicated during the verification proces, I filter for unique species names.
#   # verified_species_frame <- verified_species_frame[!duplicated(verified_species_frame$verifiedSpeciesName),]
# 
# # Pasting in the genus name and rewriting the hybrid variable to true/false   
#   verified_species_frame$verifiedSpeciesName <- paste("Saxifraga", verified_species_frame$verifiedSpeciesName)
#   verified_species_frame$speciesName <- paste("Saxifraga", verified_species_frame$speciesName)
#   verified_species_frame$hybrid <- ifelse(verified_species_frame$hybrid == "×", TRUE, FALSE)
#   
#   
#   
# # Looking up remaining species in dat and making new variables containing the new names, associated plant_name_IDs, and hybrid information 
#   dat$verifiedSpeciesName <- NA
#   dat$plantNameID <- NA
#   dat$hybrid <- FALSE
#   dat$verifiedGenus <- NA
#  
# # Appending new information to dat   
#   for (i in 1:nrow(verified_species_frame)){
#     dat$verifiedSpeciesName[dat$species == verified_species_frame[i,3]] <- verified_species_frame[i,1]
#     dat$plantNameID[dat$species == verified_species_frame[i,3]] <- verified_species_frame[i,2]
#     dat$hybrid[dat$species == verified_species_frame[i,3]] <- verified_species_frame[i,4]
#     dat$verifiedGenus[dat$species == verified_species_frame[i,3]] <- verified_species_frame[i,5]
#   }
#   
# Verified species names, their associated plant IDs, and hybrid information have now been appended to dat, and we can proceed.

# Filtering out all NA plant_name_IDs 

dat <- dat %>%
  filter(
    !is.na(plant_name_ID)
  )  




# The digital elevation model is provided in the biogeo library. Its resolution is 10 minutes, which is about 20kmx20km per pixel.
dem <- raster("wc2.1_2.5m_elev//wc2.1_2.5m_elev.tif", xmn=-180, xmx=180, ymn=-60, ymx=90)

# Moving points that fall outside the dem raster  

dat_biogeo <- keepmainfields(dat, 
                             ID = "ID", 
                             Species = "species", 
                             x = "decimallongitude", 
                             y = "decimallatitude"
)


dat_biogeo <- nearestcell(dat_biogeo, dem)[[1]]
# test <- test[[1]]

# Making a clone of dat and turning it to sf  
dat_sf <- dat_biogeo %>%
  st_as_sf(coords = c("x", "y")
  ) %>%
  as_Spatial()
crs(dat_sf) <- crs(dem)


elevation_data <- raster::extract(dem, dat_sf, method = 'simple', df = TRUE)
names(elevation_data) <- c("ID", "elevationData")

elevation_data <- elevation_data[-1]

# elev_points <- vector()
# 
# # Running script to look up elevation values in the dem.
# for (i in 1:nrow(dat)){
#   elev_extract <- raster::extract(dem, dat[4:5][i,], method = 'simple', df = TRUE)[[2]]
#   elev_points[i] <- elev_extract
#   print(paste(i, dat$iso3[i]))
# }
# 
dat <- cbind(dat, elevation_data)


dat <- dat %>%
  transmute(ID = seq.int(nrow(dat)),
            speciesName = species,
            decimallongitude,
            decimallatitude,
            coordinateUncertaintyInMeters =  as.numeric(coordinateUncertaintyInMeters),
            elevationData,
            year,
            basisOfRecord,
            countryCode,
            taxonRank,
            gbifID,
            plantNameID = plant_name_ID,
            genus
  )


# saving data

# date: 19/10
write.csv(dat, file = "data//raw//all_saxifrage_occurrences_19_10.csv", row.names = FALSE)


