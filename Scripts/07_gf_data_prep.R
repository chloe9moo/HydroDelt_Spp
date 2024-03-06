## PREP OCCURRENCE DATA + ENVIRONMENTAL DATA FOR MODELS

library(tidyverse); library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

#read in site data
# file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
# occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
#                                                          long = col_number(),
#                                                          site_id = col_character(),
#                                                          COMID = col_character(),
#                                                          gage_no_15yr = col_character(),
#                                                          dist2gage_m_15yr = col_number(),
#                                                          dist2strm_m_flw = col_number()))

g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv")) %>%
  select(-c(station_nm, site_tp_cd, coord_acy_cd, dec_coord_datum_cd, alt_va, alt_acy_va, alt_datum_cd, data_type_cd,
            parm_cd, stat_cd, ts_id, count_nu, period, int_high))

#prep env dat ----
#load in hit metrics
hit <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240130.csv"))

##determine sites needing adjusted HIT metrics ----
#intermittent streams and streams with zero flow days
adj.sites <- data.frame(site_no = hit$site_no, is_intermit = rep(NA, nrow(hit)), is_zero_flw = rep(NA, nrow(hit)), prop_zero_flw = rep(NA, nrow(hit)))

for(i in seq_len(nrow(adj.sites))) {
  site <- adj.sites$site_no[i]
  
  #check if flow type is intermittent
  adj.sites[i, ]$is_intermit <- g.info[g.info$site_no == site, ]$flw_type == "Int"
  
  #check for zero flow days
  #original raw flow dat
  f <- file.path(PATH, "02_EnvDat/raw_daily_flow", 
                 paste0(site, ifelse(file.exists(file.path(PATH, "/02_EnvDat/raw_daily_flow", paste0(site, "_long.csv"))), #check if clipped file exists
                                     "_long.csv", ".csv")))
  q <- read_csv(f, col_types = cols(site_no = col_character(), Date = col_date(), X_00060_00003 = col_number())) %>%
    mutate(month = format(Date, "%m"), day = format(Date, "%d"))
  
  if(!any(names(q) %in% "X_00060_00003")) { next } #skip those that didn't have the right flow type
  
  adj.sites[i, ]$is_zero_flw <- any(q$X_00060_00003 == 0, na.rm = TRUE) #any zero flw days?
  
  adj.sites[i, ]$prop_zero_flw <- round((sum(q$X_00060_00003 == 0, na.rm = TRUE) / nrow(q)) * 100, 3) #how many zero flow days?
  
  rm(f, q)
}

rm(i, site)

#set up hit metrics, keep adjusted for appropriate sites in df
adj_vars <- gsub("_adj", "", names(hit)[grepl("_adj", names(hit))])
adj.sites <- filter(adj.sites, is_intermit == TRUE | is_zero_flw == TRUE)

new.hit <- hit[hit$site_no %in% adj.sites$site_no, ] %>%
  select(-all_of(adj_vars)) %>%
  rename_with(~gsub("_adj", "", .x))
  
hit <- bind_rows(
  hit[!hit$site_no %in% new.hit$site_no, names(hit)[!grepl("_adj", names(hit))]],
  new.hit
)

#bind to info
env.dat <- left_join(g.info, hit)

rm(new.hit, adj.sites, adj_vars)

## add hydro alt (McManamay et al. 2022) ----
h.alt <- read_csv(paste0(PATH, "/02_EnvDat/hydro_alt_disturb/predicted_alteration_US_model.csv"))
h.alt <- h.alt %>% select(COMID, contains("pn"))

env.dat <- left_join(env.dat, h.alt) 

#fix NAs (no COMID in hydro alt dataset)
env.dat <- env.dat %>% mutate(index = row_number())

na.h.alt <- filter(env.dat, is.na(pnSeasonal)) %>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269, remove = FALSE) %>%
  st_transform(5070)
na.coms <- na.h.alt$COMID

nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% #nhd streams prev. clipped to region (good to use clipped because of MEM usage)
  st_transform(5070)
  
# Use a while loop to repeat na replacement until no NAs (search for nearest stream with h.alt data)
while (nrow(na.h.alt) != 0) {
  message(nrow(na.h.alt), " streams with no matching alteration stream ID...")
  
  #redo stream snap to find next nearest stream
  nhd.na <- filter(nhd, !COMID %in% na.coms) #filter out NA comids
  
  nhd_near <- st_nearest_feature(na.h.alt, nhd.na, check_crs = TRUE) #identify nearest nhd stream
    
  nhd_near <- nhd.na[nhd_near, ] #pull sites
    
  na.h.alt$COMID_atlNA <- nhd_near$COMID #add COMID to data
    
  dist <- as.vector(st_distance(st_transform(na.h.alt, 5070), nhd_near, by_element = TRUE)) #get distance to nearest nhd strm + add to dataset
  na.h.alt$dist2strm_m_nhd_NA <- round(dist, 3)
    
  #rejoin w/ alt data and find if still NAs 
  tmp.na <- na.h.alt %>% 
    select(!contains(names(h.alt)), contains("COMID")) %>% 
    st_drop_geometry() %>%
    left_join(h.alt, by = c("COMID_atlNA" = "COMID"))
  env.dat <- bind_rows(filter(env.dat, !index %in% tmp.na$index), tmp.na)
    
  #check for NAs still
  na.h.alt <- filter(env.dat, is.na(pnSeasonal)) %>%
    st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269, remove = FALSE) %>%
    st_transform(5070)
  na.coms <- c(na.coms, na.h.alt$COMID_atlNA)
  
  }
  
env.dat <- env.dat %>% 
  # select(-index) %>%
  select(-c(COMID_atlNA, dist2strm_m_nhd_NA)) #good to check distance before removing!!

rm(nhd.na, na.coms, na.h.alt, tmp.na, dist, nhd_near, nhd)

##add stream temp ----
strm_temp <- read_csv(paste0(PATH, "/02_EnvDat/predicted_stream_temps_summ_hist.csv"))

env.dat <- left_join(env.dat, strm_temp)

write_csv(env.dat, paste0(PATH, "/02_EnvDat/env_dat_combined_allsites.csv"))

# some plots ----
