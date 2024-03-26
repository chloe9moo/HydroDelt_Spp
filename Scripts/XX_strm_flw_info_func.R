## FUNCTION TO ATTACH FLOW CLASS + NHD COMID TO DATA ##

if(!"sf" %in% (.packages())) { #needs sf package to work
  message("package sf not loaded, loading now...")
  library(sf)
}

if(!"tidyverse" %in% (.packages())) { #needs sf package to work
  message("package sf not loaded, loading now...")
  library(tidyverse)
}

get_strm_flow_info <- function(site_data, 
                               id_col = c("site_no"), 
                               coord_cols = c("Long", "Lat"),
                               site_crs = 4269,
                               get_COMID = TRUE,
                               get_COMID_dist = TRUE,
                               get_flw_type = TRUE,
                               get_flw_dist = TRUE) {
  
  #make sites spatial
  tmp <- site_data %>% select(matches(paste0(c(id_col, coord_cols), collapse = "|"))) %>% distinct()
  tmp <- tmp[!apply(is.na(tmp[, coord_cols]), 1, all), ] #remove sites with no spatial info
  
  sf.sites <- st_as_sf(tmp, coords = coord_cols, crs = site_crs)
  sf.sites <- st_transform(sf.sites, crs = 5070)
  
  if(get_COMID) {
    
    #bring in NHD streams, transform to projected CRS for getting measurements
    nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp"))
    nhd <- st_transform(nhd, crs = 5070) #conus albers NAD83
    
    #assign COMID to each site
    nr.line <- st_nearest_feature(sf.sites, nhd, check_crs = TRUE) #find nearest nhd strm
    nhd_near <- nhd[nr.line,]
    tmp$COMID <- nhd_near$COMID #add comid match to sites
    
    #get distance to the nearest COMID
    if(get_COMID_dist) {
      dist <- as.vector(st_distance(sf.sites, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm
      tmp$dist2strm_m_nhd <- round(dist, digits = 2)
      
    }
    
    message("COMID for sites retrieved...")
    
  }

  if(get_flw_type) {
    
    #read in classified streams
    flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% 
      select(-PopupInfo) %>%
      st_transform(5070)
    flw <- flw %>% filter(Flow_type != "BigR")
    
    #assign flow to nearest site
    nr.line <- st_nearest_feature(sf.sites, flw, check_crs = TRUE) #find nearest strm not river
    flw_near <- flw[nr.line,]
    tmp$flw_name <- flw_near$Name #add strm names
    tmp$flw_type <- flw_near$Flow_type #add flow types
    
    if(get_flw_dist) {
      
      dist <- as.vector(st_distance(sf.sites, flw_near, by_element = TRUE)) #get distance from site to strm assigned
      tmp$dist2strm_m_flw <- round(dist, digits = 2)
      
    }
    
    message("Flow class for sites retrieved...")

  }
 
  site_data <- left_join(site_data, tmp)
  
  return(site_data)
  
}