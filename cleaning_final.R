
# libraries
  library(tidyverse) ## bibliotek med mange gode funktioner
  library(CoordinateCleaner) ## bibliotek til rensefunktioner
  library(raster)
  library(biogeo)
  library(dismo)
  library(sf)


# Creating climate database, I need this to look up observations without associated climate data.

  climate <- getData(name = "worldclim", 
                   var = "bio", 
                   res = 2.5)


# Inloading distribution matrix, and species list, which are necessary for cleaning coordinates based on known species ranges.
  comm <- readRDS("data//Shapefiles//comm.rds")

# Inloading shapefile for the world. This is needed for range filtering.
  world <- st_read("data//Shapefiles//level3.shp")




##### 
#
#  Loading raw presence data. 
#
#####

  dat_raw <- read.csv("data//raw//all_saxifrage_occurrences_19_10.csv")

# Selecting variables, and doing an initial filtering. 
# The filtering is to remove any data without coordinate values, any hybrid species, and data without the desired type of record.     
  dat_raw <- dat_raw %>% 
    transmute(ID,
            speciesName, 
            decimallongitude, 
            decimallatitude,
            coordinateUncertaintyInMeters,
            elevationData,
            year,
            basisOfRecord, 
            countryCode,
            taxonRank,
            genus,
            plantNameID,
            gbifID
  ) %>% 
  filter(!is.na(decimallongitude), 
         !is.na(decimallatitude),
         !is.na(plantNameID),
         basisOfRecord == "MACHINE_OBSERVATION" |
           basisOfRecord == "HUMAN_OBSERVATION" |
           basisOfRecord == "OBSERVATION" |
           basisOfRecord == "PRESERVED_SPECIMEN" |
           basisOfRecord == "LIVING_SPECIMEN",
         # hybrid == FALSE,
         plantNameID %in% colnames(comm)
  )        

# Filtering out duplicates     
  dat_raw <- cc_dupl(dat_raw, lon = "decimallongitude", lat = "decimallatitude", species = "speciesName", value = "clean")



########## Range cleaning ####
# This part filters out all points that fall outside the range shapefiles. 

  id_list <- as.data.frame(cbind(dat_raw$speciesName, dat_raw$plantNameID))
  names(id_list) <- c("speciesName", "plantNameID")

  id_list <- id_list[!duplicated(id_list$speciesName),]
  id_list <- id_list[id_list$plantNameID %in% colnames(comm),]

  dat_raw_intersected <- dat_raw[FALSE, ]

for (i in 1:nrow(id_list)){
  # if(is.na(id_list_test[i,2]) == TRUE){next}
  print(id_list[i,1])
  country_list <- comm[, id_list[i,2]]
  
  # if(exists(paste0(country_list))){next}
  
  country_names <- names(country_list[country_list == 1])
  
  world_subset <- world %>%
    filter(LEVEL_3_CO %in% country_names)
  
  dat_subset <- dat_raw[dat_raw$speciesName == id_list[i,1], ]
  subset_matrix <- matrix(c(dat_subset$decimallongitude, dat_subset$decimallatitude), ncol = 2)
  
  subset_multi = st_multipoint(x = subset_matrix)
  subset_p = st_cast(st_sfc(subset_multi), "POINT")
  st_crs(subset_p) <- st_crs(world_subset)
  
  inter <- st_intersects(subset_p, world_subset, sparse = FALSE)
  iter <- vector()
  
  for (j in 1:nrow(inter)){
    iter[j] <- ifelse(TRUE %in% inter[j,], TRUE, FALSE )
  }
  
  dat_raw_intersected <- rbind(dat_raw_intersected, dat_subset[iter,])
  rm(country_list)
  
}

dat_raw <- dat_raw_intersected   
#     
##### End of range cleaning    

###### Manual cleaning
#Adding columns needed for the biogeo library functions
  dat_raw <- keepmainfields(dat_raw,
                          ID = "ID",
                          Species = "speciesName",
                          x = "decimallongitude",
                          y = "decimallatitude",
                          others = c("coordinateUncertaintyInMeters",
                                     "elevationData",
                                     "year",
                                     "basisOfRecord",
                                     "countryCode",
                                     # "hybrid",
                                     "taxonRank",
                                     "genus",
                                     "plantNameID",
                                     "gbifID"
                          )
)

#  Points at coastlines may fall outside the climate raster, due to its resolution. 
#  These points can be moved, if a neighbour rastercell has a value.    
dat_raw <- nearestcell(dat_raw, climate[[1]])[[1]]

# Flagging all points without climate data that were not moved.
dat_raw <- missingvalsexclude(climate[[1]], dat_raw)

# Automatic cleaning using CoordinateCleaner library
flags_cc <- clean_coordinates(x = dat_raw, lon = "x", lat = "y", species = "Species",
                              tests = c(
                                "capitals",
                                "centroids",
                                "equal",
                                "gbif",
                                "institutions",
                                "zeros"
                              ),
                              value = "flagged"
)


# Appending cc flags to the dataframe    
dat_raw <- cbind(dat_raw, flags_cc)


# Removing values flaggd by either biogeo or coordinate cleaner. Biogeo sets all bad points to 1 in the exclude column,
# while coordinate cleaner sets all bad points to TRUE in the flags_cc column.    
dat_raw <- dat_raw %>%
  filter(Exclude == 0,
         flags_cc == TRUE
  )

# Removing unneeded columns 
dat_raw <- dat_raw %>%
  transmute(Species,
            decimallongitude = x,
            decimallatitude = y,
            gbifID,
            coordinateUncertaintyInMeters,
            year,
            basisOfRecord,
            countryCode,
            elevationData,
            taxonRank,
            plantNameID
  )

#Saving file 
write.csv(dat_raw, file = "data//Species//all_saxifrage_occurrences_cleaned_19_10.csv", row.names = FALSE)

