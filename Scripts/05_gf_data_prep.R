## PREP OCCURRENCE DATA + ENVIRONMENTAL DATA FOR MODELS

library(tidyverse); library(sf)

PATH <- getwd()

#bring in occurrence data (for not this will be the already filtered datasets)
#read in occ data ----
fish <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/ARMOOK_Fishes_bysite23a_StreamCat_Flow_Full_15km.csv")) %>% select(-`...1`) %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 26915, remove = FALSE) #NAD83 UTM Zone 15
bug <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/BenthicInsect_23_StreamCat_Flow_Full_15km.csv")) %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 4269, remove = FALSE) #NAD83 LL

b.tax <- list(fish, bug)

# snap to nearest nhd stream ----
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% 
  st_transform(5070) #albers, projected crs

b.tax <- lapply(b.tax, function(x) { 
  
  nhd_near <- x %>%
    st_transform(5070) %>% #NAD83 Conus Albers
    st_nearest_feature(., nhd, check_crs = TRUE) #identify nearest nhd stream
  
  nhd_near <- nhd[nhd_near, ] #pull sites
  
  x$COMID_vCEM <- nhd_near$COMID #add COMID to data (vCEM to keep separate from original columns, update once cleaned data is used)

  dist <- as.vector(st_distance(st_transform(x, 5070), nhd_near, by_element = TRUE)) #get distance to nearest nhd strm + add to dataset
  x$dist2strm_m_nhd <- round(dist, 3)
  
  return(x)
  
})

#note! luckily seems to match original comid for the most part (like 10 differences)

# add hydro alt (McManamay et al. 2022) ----
h.alt <- read_csv(paste0(PATH, "/02_EnvDat/hydro_alt_disturb/predicted_alteration_US_model.csv"))
h.alt <- h.alt %>% select(COMID, contains("pn"))

lapply(b.tax, function(x) {
  
  tmp <- x %>%
    left_join(h.alt, by = c("COMID_vCEM" = "COMID"))
  
  #fix NAs (no COMID in hydro alt dataset)
  na <- tmp %>% filter(is.na(pnHA_rank))
  nhd.na <- filter(nhd, !COMID %in% na$COMID_vCEM)
  
  #redo stream snap to find next nearest stream
  nhd_near <- na %>%
    st_transform(5070) %>% #NAD83 Conus Albers
    st_nearest_feature(., nhd.na, check_crs = TRUE) #identify nearest nhd stream
  
  nhd_near <- nhd.na[nhd_near, ] #pull sites
  
  na$COMID_vCEM_NA <- nhd_near$COMID #add COMID to data (vCEM to keep separate from original columns, update once cleaned data is used)
  
  dist <- as.vector(st_distance(st_transform(na, 5070), nhd_near, by_element = TRUE)) #get distance to nearest nhd strm + add to dataset
  na$dist2strm_m_nhd_NA <- round(dist, 3)
  
  #rejoin 
  tmp.na <- na %>% 
    select(!contains(names(h.alt)), contains("COMID")) %>% 
    left_join(h.alt, by = c("COMID_vCEM_NA" = "COMID"))
  
  #TO DO NEXT: FIGURE OUT WHETHER/HOW TO LOOP THROUGH NAS AND REPLACE WITH NEAREST NHD STRM W/ HYDRO ALT DATA
  
})




library(mapview)
mapviewOptions(fgb = FALSE)

mapview(nhd.na) + mapview(na)
#add stream temp (cv + mean monthly + seasonal monthly)














