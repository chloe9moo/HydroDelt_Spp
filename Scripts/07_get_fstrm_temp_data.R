## OBTAIN AND PREP DATA FOR FINE SCALE STREAM TEMP

## STEPS:
#1. Download temperature data for USGS gages
#2. Calculate monthly mean temperature for each month of data (weekly?)
#3. Compile external water temp sources
#4. Exclude months with fewer than 20 days temp data
#5. Collect monthly (weekly?) mean air temp from PRISM
#6. Extract air temp at sites

#OLD:
#4. Collect monthly mean air temp from the National Oceanic and Atmospheric Agency National Center for Climate Data
##4a. From the county where each site was located for a corresponding time period to the water temp data

library(sf); library(tidyverse); library(maps)
library(dataRetrieval)#; library(httr2)
options(readr.show_col_types = FALSE)

PATH <- getwd()

# source(paste0(PATH, "/Scripts/XX_api_token.R"))

#OLD
# hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv")) #<< all stations in ARMOOK, HIT calculated
# g.info <- dataRetrieval::whatNWISdata(siteNumbers = unique(hit$STAID_0)) #get xy for all sites and save for usgs sites for later
# write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_all_prevHITsites_info.csv"))
#NEW
# hit <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"))

#1. Download temperature data for ARMOOK USGS gages ----
# https://waterdata.usgs.gov/blog/dataretrieval/

#00010 = pCode for daily temp in celsius
temp.gage <- data.frame()
for(state_nm in c("AR", "MO", "OK")) {
  tmp <- whatNWISdata(stateCd = state_nm,
                      parameterCd = "00010",
                      service = "dv",
                      statCd = "00003")
  temp.gage <- rbind(temp.gage, tmp)
}

temp.gage <- temp.gage %>%
  mutate(period = as.Date(end_date) - as.Date(begin_date))

#save for later (JIC)
write_csv(temp.gage, paste0(PATH, "/02_EnvDat/raw_stream_temp/usgs_gage_wtemp_ARMOOK_info.csv"))
temp.gage <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/usgs_gage_wtemp_ARMOOK_info.csv"))

temp.gage <- temp.gage  %>%
  #remove lake/spring data
  filter(!site_tp_cd %in% c("LK", "SP")) %>%
  #at least one year of data, more than 20 records (exclude months with less than 20 days of data)
  filter(period >= 365 & count_nu >= 20) %>%
  #not max. or min. stat data
  filter(is.na(stat_cd) | stat_cd == "00003")

#get daily values temp data
temp.dv <- readNWISdv(siteNumbers = unique(temp.gage$site_no), 
                      parameterCd = "00010",
                      startDate = min(temp.gage$begin_date),
                      endDate = max(temp.gage$end_date))

write_csv(temp.dv, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_", gsub("-", "", Sys.Date()), ".csv"))

#deal with multi-location across stream record sites
temp.dv <- temp.dv %>% mutate(id_no = row_number()) #for checking row loss later

tmp <- temp.dv[!is.na(temp.dv$X_.Left.Bank._00010_00003) | !is.na(temp.dv$X_.Right.Bank._00010_00003) | !is.na(temp.dv$X_Left.Bank_00010_00003), ]
tmp <- tmp %>% mutate(X_00010_00003 = rowMeans(select(., X_.Left.Bank._00010_00003, X_.Right.Bank._00010_00003, X_Left.Bank_00010_00003), na.rm = T))

tmp1 <- temp.dv %>% 
  filter(if_any(contains("Cross.section"), ~!is.na(.x))) %>%
  mutate(X_00010_00003 = rowMeans(select(., ends_with("_00003"), -X_00010_00003), na.rm = T),
         X_00010_00003 = round(X_00010_00003, digits = 1))

tmp2 <- temp.dv %>%
  filter(!id_no %in% c(tmp$id_no, tmp1$id_no))

tmp3 <- bind_rows(tmp, tmp1, tmp2)

# View(filter(temp.dv, !id_no %in% tmp3$id_no))
# View(filter(tmp3, is.na(X_00010_00003)))

temp.dv <- tmp3 %>%
  rename(w_temp = X_00010_00003) %>%
  select(agency_cd, site_no, Date, w_temp) %>%
  mutate(date = as.Date(Date),
         month = month(Date),
         day = day(Date),
         year = year(Date),
         week = week(date))

write_csv(temp.dv, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_", gsub("-", "", Sys.Date()), ".csv"))

rm(list = ls(pattern = "tmp"))

#2. Calculate weekly+monthly mean temperature for each month of data ----
##summarizing function
get_month_week_temp_sum <- function(daily_temp_dat, temp_col = c("w_temp", "mn_daily_temp")) {
  
  return_list <- vector("list", 2)
  
  t_name <- temp_col[temp_col %in% names(daily_temp_dat)]
  
  #MONTHLY
  return_list[[1]] <- daily_temp_dat %>%
    group_by(site_no, year, month) %>%
    summarise(mn_monthly_temp = mean(get(t_name), na.rm = T),
              ct = n()) %>%
    arrange(year, month, site_no) %>%
    ungroup()
  
  #WEEKLY
  return_list[[2]] <- daily_temp_dat %>%
    group_by(site_no, year, week) %>%
    summarise(mn_weekly_temp = mean(get(t_name), na.rm = T),
              ct = n()) %>%
    arrange(year, week, site_no) %>%
    ungroup()
  
  names(return_list) <- c("monthly", "weekly")
  
  return(return_list)
}

#load daily data
temp.all <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_20240312.csv"))

temp.all <- temp.all %>%
  filter(!is.na(w_temp) & w_temp < 60) #remove missing data or erroneous data (see Hare et al. 2023)

#summarize
summ.temp <- get_month_week_temp_sum(temp.all)

summ.temp <- lapply(summ.temp, function(x) {
  x %>% mutate(across(contains("mn_"), ~ round(.x, digits = 2))) #match precision of original data
})

write_csv(summ.temp[["monthly"]], paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_monthly_temp_ARMOOK_all_", gsub("-", "", Sys.Date()), ".csv"))
write_csv(summ.temp[["weekly"]], paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_weekly_temp_ARMOOK_all_", gsub("-", "", Sys.Date()), ".csv"))

rm(temp.all, summ.temp, hit)

#3+. Clean, compile, collected temperature data ----
dirs <- list.dirs(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data"), recursive = FALSE)

tmp <- list.files(dirs[[1]], "[0-9].csv", full.names = TRUE)

##3a. load in external temp sources ----
###eel ----
h <- readxl::read_excel(list.files(path = dirs[[1]], pattern = "[Hh]eader", full.names = TRUE), skip = 6) #get site info
t <- lapply(list.files(dirs[[1]], "[0-9].csv", full.names = TRUE), read_csv, skip = 1) #read in individual temp files
#one was wonky, so I corrected the column names in command line and reload a different way
t[[8]] <- read_delim(list.files(dirs[[1]], "[0-9].csv", full.names = TRUE)[[8]], delim = "\t", escape_double = FALSE) 

for(i in 1:length(t)){
  colnames(t[[i]]) <- c("row_id", "date", "temp") #rename column for consistency
  t[[i]]$site_no <- gsub(".csv", "", list.files(dirs[[1]], ".csv")[[i]])
  
  t[[i]] <- t[[i]] %>% 
    mutate(date = sub(" .*$", "", date),
           n.c = nchar(date))
  if(all(unique(t[[i]]$n.c) == 8)) { #some files have diff date format
    t[[i]] <- t[[i]] %>%
      mutate(date = as.Date(date, format = "%m/%d/%y"),
             year = year(date),
             month = month(date),
             week = week(date),
             day = day(date)) %>%
      select(-n.c)
  } else {
    t[[i]] <- t[[i]] %>%
      mutate(date = as.Date(date, format = "%m/%d/%Y"),
             year = year(date),
             month = month(date),
             week = week(date),
             day = day(date)) %>%
      select(-n.c)
  }
}

t <- bind_rows(t)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Eel_Header__Data/eel_water_temp_all_comped.csv"))
# t <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Eel_Header__Data/eel_water_temp_all_comped.csv"))

t <- t %>% 
  select(-row_id) %>%
  group_by(site_no, year, month, week, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% 
  ungroup() %>%
  left_join(., h %>% select(FileName, Lat, Long), by = c("site_no" = "FileName")) %>%
  select(-site_no) %>% #noticed that there is overlap in sites, so going to just rename each site
  group_by(Lat, Long) %>%
  mutate(site_no = paste0("eel_", cur_group_id())) %>%
  arrange(site_no, year, month, day)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/eel_water_daily_temp.csv"))

###MO JODI----
h <- read_csv(list.files(path = dirs[[2]], pattern = "[Hh]eader", full.names = TRUE))
t <- lapply(list.files(dirs[[2]], "[0-9].csv", full.names = TRUE), read_csv) #read in individual temp files

t <- lapply(t, function(x) {
  df <- x %>% select(SiteID, Date2, Temp_C)
  colnames(df) <- c("site_no", "date", "temp")
  df <- df %>% 
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           year = year(date),
           month = month(date),
           week = week(date),
           day = day(date)) %>%
    group_by(site_no, year, month, week, day) %>%
    summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% 
    ungroup()
  return(df)
})

t <- bind_rows(t)
t <- t %>% 
  left_join(., h %>% select(SiteID, Lat, Long), by = c("site_no" = "SiteID")) %>% #NOTE, SOME SITES DON'T HAVE LAT LONG DATA ??? FOR SOME REASON!
  select(-site_no) %>%
  group_by(Lat, Long) %>%
  mutate(site_no = paste0("mo_jodi_", cur_group_id())) %>%
  arrange(site_no, year, month, day) %>%
  distinct()
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/MO_water_daily_temp.csv"))

### strawberry ----
h <- readxl::read_excel(list.files(path = dirs[[3]], pattern = "[Hh]eader", full.names = TRUE), skip = 6)
t <- lapply(list.files(dirs[[3]], "[0-9].xlsx", full.names = TRUE), readxl::read_excel) #read in individual temp files

for(i in 1:length(t)){
  colnames(t[[i]]) <- c("date", "temp") #rename column for consistency
  t[[i]]$site_no <- gsub(".xlsx", "", list.files(dirs[[3]], "[0-9].xlsx")[[i]])
  
  t[[i]] <- t[[i]] %>%
    mutate(date = sub(" .*$", "", date),
           date = as.Date(date, format = "%Y-%m-%d"),
           year = year(date),
           month = month(date),
           week = week(date),
           day = day(date))
}

t <- bind_rows(t) 
#fix misspelling
t <- t %>%
  mutate(site_no  = ifelse(site_no == "10996600__Straw_HulettR_9.23.2019", "10996600__Straw_HulettRd_9.23.2019", site_no))
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Strawberry_Header__Data/strawb_water_temp_all_comped.csv"))
# t <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Strawberry_Header__Data/strawb_water_temp_all_comped.csv"))

#daily temp
t <- t %>% 
  group_by(site_no, year, month, week, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>%
  ungroup() %>% 
  left_join(., h %>% select(FileName, Lat, Long), by = c("site_no" = "FileName")) %>%
  select(-site_no) %>% #noticed that there is overlap in sites, so going to just rename each site
  group_by(Lat, Long) %>%
  mutate(site_no = paste0("strawb_", cur_group_id())) %>%
  arrange(site_no, year, month, day)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/strawb_water_daily_temp.csv"))

###adams lab ----
h <- read_csv(list.files(path = dirs[[4]], pattern = "[Hh]eader", full.names = TRUE), skip = 6)
h <- h %>% mutate(FileName = gsub(".csv", "", FileName))
t <- lapply(list.files(dirs[[4]], ".csv", full.names = TRUE)[-22], read_csv) #read in individual temp files, ignore header file

for(i in 1:length(t)){
  t[[i]] <- t[[i]] %>% select(Date, Temperature_C) %>% rename(date = Date, temp = Temperature_C) #rename column for consistency
  t[[i]]$site_no <- gsub(".csv", "", list.files(dirs[[4]], ".csv")[-22][[i]])
  
  t[[i]] <- t[[i]] %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           year = year(date),
           month = month(date),
           week = week(date),
           day = day(date))
}

t <- bind_rows(t)

#daily temp
t <- t %>% 
  group_by(site_no, year, month, week, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., h %>% select(FileName, Lat, Long), by = c("site_no" = "FileName")) %>%
  select(-site_no) %>% #noticed that there is overlap in sites, so going to just rename each site
  group_by(Lat, Long) %>%
  mutate(site_no = paste0("uca_", cur_group_id())) %>%
  arrange(site_no, year, month, day)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/uca_water_daily_temp.csv"))

rm(h, t, summ.temp)

##3b. compile ALL data into single file, with location info ----
#dailies
#usgs dat
us <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_20240312.csv"))
temp.gage <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/usgs_gage_wtemp_ARMOOK_info.csv"))

us <- left_join(us, 
                temp.gage %>% 
                  rename(Lat = dec_lat_va, Long = dec_long_va) %>% 
                  select(site_no, Lat, Long) %>%
                  distinct()) %>%
  select(-agency_cd, -Date, -date) %>%
  rename(mn_daily_temp = w_temp)

#non-usgs dat
nu <- lapply(list.files(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data"), "daily", full.names = TRUE), read_csv)

nu <- lapply(nu, function(x) { 
  x %>% mutate(site_no = as.character(site_no))
})

nu <- bind_rows(nu)

wtemp <- bind_rows(nu, us)
write_csv(wtemp, paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_daily_temp_", gsub("-", "", Sys.Date()), ".csv"))

##get compiled monthly + weekly summaries
all.temp <- get_month_week_temp_sum(wtemp)

#add lat long back in
site.info <- wtemp %>% select(site_no, Lat, Long) %>% distinct()
all.temp <- lapply(all.temp, function(x) {
  x %>% left_join(., site.info)
})

write_csv(all.temp[["monthly"]], paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_temp_", gsub("-", "", Sys.Date()), ".csv"))
write_csv(all.temp[["weekly"]], paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_weekly_temp_", gsub("-", "", Sys.Date()), ".csv"))

rm(all.temp, nu, us, wtemp, site.info)

##3c. Attach gage data to streams + add flow cat ----
source(paste0(PATH, "/Scripts/XX_strm_flw_info_func.R"))

temp.site.info <- get_strm_flow_info(site_data = wtemp[["monthly"]])
temp.site.info <- temp.site.info %>% select(site_no, Lat, Long, COMID, contains("dist"), contains("flw")) %>% distinct()

write_csv(temp.site.info, paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_location_info.csv"))

#4. Minimum data exclusion ----
f <- list.files(paste0(PATH, "/02_EnvDat/raw_stream_temp"), "comb_water", full.names = TRUE)
f <- f[!grepl("location_info", f)]
wtemp <- lapply(f, read_csv)
names(wtemp) <- c("daily", "monthly", "weekly")

# check <- lapply(wtemp, function(x){
#   x %>% filter(is.na(Lat))
# })

#minimum for monthly is at least 20 days of data per included month
wtemp[["monthly"]] <- wtemp[["monthly"]] %>% 
  filter(ct >= 20) %>%
  mutate(date = with(., sprintf("%d-%02d", year, month))) %>%
  filter(!is.na(Lat))

# #minimum for weekly is at least 5 days of data per included week
# wtemp[["weekly"]] <- wtemp[["weekly"]] %>%
#   filter(ct >= 5) %>%
#   filter(!is.na(Lat))

#5. Get prism data across ALL dates ----
library(prism); library(terra)

#load hit dates 
hit.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
hit.info <- hit.info[names(hit.info) %in% c("site_no", "dec_long_va", "dec_lat_va", "hit_start_date", "hit_end_date")]

#set prism download directory
prism_set_dl_dir(paste0(PATH, "/02_EnvDat/raw_air_temp/prism"))

#for now at least, going to go with monthly for space usage
# date_vect <- wtemp[["daily"]] %>% mutate(date = as.Date(paste0(year, "-", month, "-", day))) %>% arrange(date)
# date_vect <- unique(date_vect[date_vect$year > 1980, ]$date)
# get_prism_dailys(type = "tmean", dates = date_vect[2:4])

date_vect <- wtemp[["monthly"]]$year
min <- min(c(date_vect, year(hit.info$hit_start_date)))
max <- max(c(date_vect, year(hit.info$hit_start_date)))
get_prism_monthlys(type = "tmean",
                   years = min:max,
                   mon = seq(1, 12),
                   keepZip = TRUE, #can remove the unzipped files after clipping
                   keep_pre81_months = TRUE)

#get lat long of all gages + bio sites for clipping
occ.sites <- read_csv(paste0(PATH, "/01_BioDat/sites_alltax_inthigh_allinfo_20240717.csv")) %>% 
  select(long, lat) %>% rename(Long = long, Lat = lat) %>% distinct()
coords <- wtemp[["monthly"]] %>% select(Lat, Long) %>% distinct()
coords <- bind_rows(
  coords,
  hit.info %>% select(dec_lat_va, dec_long_va) %>% rename(Lat = dec_lat_va, Long = dec_long_va),
  occ.sites
)
coords <- st_as_sf(coords, coords = c("Long", "Lat"), crs = 4269)
coords <- st_bbox(coords)
#add some wiggle room
coords[c("xmax", "ymax")] <- coords[c("xmax", "ymax")] + 0.5
coords[c("xmin", "ymin")] <- coords[c("xmin", "ymin")] - 0.5

prism.list <- list.dirs(paste0(PATH, "/02_EnvDat/raw_air_temp/prism"), full.names = FALSE, recursive = FALSE)
prism.list <- prism.list[!prism.list %in% "clipped_rasters"] #remove folder for saving clips

for(i in seq_along(prism.list)) {
  
  rast_name <- prism.list[i]
  
  if(file.exists(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters/", rast_name, "_clipped.tif"))) { next }
  
  #load in raster
  prism.rast <- rast(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/", rast_name, "/", rast_name, ".bil"))
  #clip to extent of coordinates
  prism.c <- terra::crop(prism.rast, ext(coords))
  #save to new folder
  writeRaster(prism.c, filename = paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters/", rast_name, "_clipped.tif"), overwrite = FALSE) #using .bil to match
  
  message(round(i/length(prism.list)*100, digits = 2), "% complete")
  
}

#plot check
rast_name <- "PRISM_tmean_stable_4kmM3_200709_bil"
tmp <- rast(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters/", rast_name, "_clipped.tif"))
# # plot(prism.rast)
# # plot(prism.c)
plot(tmp)
# # # plot(hlnd$geom, add = TRUE)
plot(coords$geometry, add = TRUE)

#6. Extract temperature values for regression prep ----
#get lat long of for value extracting
occ.sites <- read_csv(paste0(PATH, "/01_BioDat/sites_alltax_inthigh_allinfo_20240717.csv")) %>%
  select(site_id_new, long_new, lat_new) %>% rename(site_no = site_id_new, Long = long_new, Lat = lat_new) %>% distinct()
coords <- wtemp[["monthly"]] %>% select(site_no, Lat, Long) %>% distinct()
coords <- bind_rows(
  coords,
  hit.info %>% select(site_no, dec_lat_va, dec_long_va) %>% rename(Lat = dec_lat_va, Long = dec_long_va),
  occ.sites
)
coords <- distinct(coords)

#get list of clipped rasters
prism.list <- list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters"))
prism.list <- prism.list[!grepl("_ppt_", prism.list)] #exclude precip data

#set dataframe
a_temp <- data.frame()
for(i in seq_along(prism.list)) {
  
  date <- gsub(".*?(\\d{6}).*", "\\1", prism.list[i])
  if(nchar(date) != 6) { next }
  
  #read in raster
  prism.rast <- rast(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters/", prism.list[i]))
  
  #extract temperature values at all gage sites
  a_temp1 <- extract(prism.rast, cbind(coords$Long, coords$Lat))
  names(a_temp1) <- "air_temp"
  a_temp1$year <- gsub("^(....).*", "\\1", date)
  a_temp1$month <- gsub(".*(.{2})$", "\\1", date)
  a_temp1 <- bind_cols(a_temp1, coords)
  
  #bind to full dataframe
  a_temp <- bind_rows(a_temp, a_temp1)
  
  message(round(i/length(prism.list)*100, digits = 2), "% complete")
  
}

write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_gage_bio_sites.csv"))

# rm(a_temp, a_temp1, coords, prism.rast, prism.c, max, min, rast_name, date, date_vect) #clean up env
# a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_strm_gage_sites.csv"))
