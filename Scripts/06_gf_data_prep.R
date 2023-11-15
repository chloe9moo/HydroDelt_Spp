## PREP OCCURRENCE DATA + ENVIRONMENTAL DATA FOR MODELS

library(tidyverse); library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

#bring in occurrence data (for not this will be the already filtered datasets)
#read in occ data ----
fish <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/ARMOOK_Fishes_bysite23a_StreamCat_Flow_Full_15km.csv")) %>% select(-`...1`) %>%
  select(-STAID_1, -STAID_2) %>%
  mutate(STAID = paste0("0", STAID)) %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 26915, remove = FALSE) #NAD83 UTM Zone 15
bug <- read_csv(paste0(PATH, "/01_BioDat/archive_occ_dat/BenthicInsect_23_StreamCat_Flow_Full_15km.csv")) %>%
  mutate(STAID = paste0("0", STAID)) %>%
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

b.tax <- lapply(b.tax, function(x) {
  
  tmp <- x %>%
    left_join(h.alt, by = c("COMID_vCEM" = "COMID")) %>%
    mutate(index = row_number())
  
  #fix NAs (no COMID in hydro alt dataset)
  na <- tmp %>% filter(is.na(pnHA_rank))
  na.coms <- na$COMID_vCEM
  nhd.na <- filter(nhd, !COMID %in% na.coms) #filter out NA comids
  
  # Use a while loop to repeat na replacement until no NAs (search for nearest stream with h.alt data)
  while (nrow(na) != 0) {
    
    cat(nrow(na), "rows in dataset\n")
    
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
    tmp <- bind_rows(filter(tmp, !index %in% tmp.na$index), tmp.na)
    
    #check for NAs still
    na <- tmp %>% filter(is.na(pnHA_rank))
    na.coms <- c(na.coms, na$COMID_vCEM_NA)
    nhd.na <- filter(nhd, !COMID %in% na.coms) #filter out NA comids

  }
  
  tmp <- tmp %>% select(-index) #remove the NA index
  
  return(tmp)
  
})


# summarize stream temp values ----
#occ data goes from 1900 to 2022, but the NCCV categorizes 1950 - 2005 as historical and 2006 - 2099 as 21st century period
temp <- read_csv(paste0(PATH, "/02_EnvDat/predicted_stream_temps_monthly.csv"))

temp.list <- list(temp %>% filter(date < "2006-01-15"), temp %>% filter(date > "2005-12-15")) #split historical and predicted data

summ.list <- lapply(temp.list, function(temp) {
  
  tmp <- temp %>%
    #prep df
    mutate(date = as.Date(date),
           year = year(date),
           month = month(date),
           season = case_when(month %in% c(12, 1, 2) ~ "winter",
                              month %in% c(3, 4, 5) ~ "spring",
                              month %in% c(6, 7, 8) ~ "summer",
                              month %in% c(9, 10, 11) ~ "fall")) %>%
    #calc summ stats for each month
    group_by(site_no, month) %>%
    mutate(across(contains("pred_strm"), list(
      mn = ~mean(.x, na.rm = T),
      med = ~median(.x, na.rm = T),
      max = ~max(.x, na.rm = T),
      min = ~min(.x, na.rm = T)
    ), .names = "mnth_{.fn}_{.col}")) %>%
    ungroup() %>%
    rename_with(.cols = contains("mnth"), ~ sub("pred_strm_temp", "temp", .x)) %>%
    #calc summ stats for each season
    group_by(site_no, season) %>%
    mutate(across(contains("pred_strm"), list(
      mn = ~mean(.x, na.rm = T),
      med = ~median(.x, na.rm = T),
      max = ~max(.x, na.rm = T),
      min = ~min(.x, na.rm = T)
    ), .names = "ssn_{.fn}_{.col}")) %>%
    ungroup() %>%
    rename_with(.cols = contains("ssn"), ~ sub("pred_strm_temp", "temp", .x)) %>%
    #get avg annual CV
    group_by(site_no, year) %>%
    mutate(across(contains("pred_strm"), ~sd(.x, na.rm = T) / mean(.x, na.rm = T), .names = "yr_cv_{.col}")) %>%
    ungroup()
  
  mnth <- tmp %>%
    select(site_no, month, contains("mnth")) %>%
    distinct() %>%
    pivot_wider(id_cols = site_no, 
                names_from = month,
                values_from = contains("mnth"),
                names_glue = "{.value}_{month}")
  ssn <- tmp %>%
    select(site_no, season, contains("ssn")) %>%
    distinct() %>%
    pivot_wider(id_cols = site_no, 
                names_from = season,
                values_from = contains("ssn"),
                names_glue = "{.value}_{season}") 
  cv <- tmp %>%
    select(site_no, contains("yr")) %>%
    distinct() %>%
    group_by(site_no) %>%
    summarize(across(contains("yr_cv"), ~ mean(.x, na.rm = T))) %>%
    rename_with(.cols = contains("yr_cv"), ~ sub("yr_cv_pred_strm_temp", "ann_temp_cv", .x))
  
  summ.temp <- left_join(mnth, ssn) %>% left_join(., cv)
  
  return(summ.temp)
})

write_csv(summ.list[[1]], paste0(PATH, "/02_EnvDat/predicted_stream_temps_summ_hist.csv"))
write_csv(summ.list[[2]], paste0(PATH, "/02_EnvDat/predicted_stream_temps_summ_future.csv"))

b.tax <- lapply(b.tax, function(df) {
  df <- df %>%
    left_join(., summ.list[[1]], by = c("STAID" = "site_no"))
})

write_csv(b.tax[[1]], paste0(PATH, "/01_BioDat/archive_occ_dat/fish_bio_env_updated_20231110.csv"))
write_csv(b.tax[[2]], paste0(PATH, "/01_BioDat/archive_occ_dat/bug_bio_env_updated_20231110.csv"))

# some plots ----

# temp %>% 
#   filter(site_no == "06906800") %>%
#   mutate(date = as.Date(date),
#          year = year(date),
#          month = as.factor(month(date))) %>%
#   filter(year < 2025 & year > 2020) %>%
#   ggplot() +
#   geom_line(aes(x = date, y = pred_strm_temp_8.5)) +
#   geom_point(aes(x = date, y = pred_strm_temp_8.5, color = month))

# g.xy <- read_csv(paste0(PATH, "/02_EnvDat/all_hit_usgs_gage_info.csv")) %>%
#   select(site_no, dec_long_va, dec_lat_va) %>% distinct() %>%
#   st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269)
# hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
#   st_transform(., crs = 4269)
# 
# tmp <- g.xy %>% left_join(., summ.list[[2]])
# 
# ggplot() +
#   geom_sf(data = hlnd) +
#   geom_sf(data = filter(tmp, !is.na(ann_temp_cv_4.5)), aes(fill = ann_temp_cv_4.5), shape = 21, size = 2.5, alpha = 0.7) +
#   geom_sf(data = filter(tmp, is.na(ann_temp_cv_4.5)), fill = "grey50", shape = 21, size = 2, alpha = 0.5) +
#   geom_sf(data = g.xy %>% filter(site_no == "07364150"), shape = 13, size = 6, color = "red") +
#   scale_fill_viridis_c(breaks = seq(0, 0.7, 0.1), limits = c(-0.001, 0.71)) +
#   guides(fill = guide_colorbar(title = "avg.\nannual\nCV")) +
#   theme_minimal()
# ggsave(paste0(PATH, "/99_figures/pred_strm_temp_cv_future_map.png"), width = 7, height = 5, bg = "white")







