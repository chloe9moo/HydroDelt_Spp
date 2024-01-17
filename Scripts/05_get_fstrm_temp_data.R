## OBTAIN AND PREP DATA FOR FINE SCALE STREAM TEMP
# 
# following the methods of Middaugh et al. 2016 (Ecology of Freshwater Fish)

## STEPS:
#1. Download temperature data for USGS gages
#2. Calculate monthly mean temperature for each month of data
#3. Exclude months with fewer than 20 days temp data
#4. Collect monthly mean air temp from the National Oceanic and Atmospheric Agency National Center for Climate Data
##4a. From the county where each site was located for a corresponding time period to the water temp data

library(sf); library(tidyverse); library(maps)
library(dataRetrieval); library(httr2)
options(readr.show_col_types = FALSE)

PATH <- getwd()

source(paste0(PATH, "/Scripts/XX_api_token.R"))

hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv")) #<< all stations in ARMOOK, HIT calculated

# g.info <- dataRetrieval::whatNWISdata(siteNumbers = unique(hit$STAID_0)) #get xy for all sites and save for usgs sites for later
# write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_all_prevHITsites_info.csv"))

#1. Download temperature data for USGS gages ----
# https://waterdata.usgs.gov/blog/dataretrieval/
#get site numbers from data:
g.nums <- unique(hit$STAID_0)

#00010 = pCode for daily temp in celsius
g.info <- whatNWISdata(siteNumbers = g.nums, parameterCd = "00010") #get gauge info
g.get <- g.info %>%
  mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
  #at least one year of data, more than 20 records (exclude months with less than 20 days of data)
  filter(period >= 365 & count_nu >= 20) %>%
  #not max. or min. stat data
  filter(is.na(stat_cd) | stat_cd == "00003")

unique(g.get$data_type_cd)

temp.dv <- readNWISdv(siteNumbers = unique(g.get[g.get$data_type_cd == "dv", ]$site_no), 
                      parameterCd = "00010",
                      startDate = min(g.get[g.get$data_type_cd == "dv", ]$begin_date),
                      endDate = max(g.get[g.get$data_type_cd == "dv", ]$end_date))
temp.dv <- temp.dv %>%
  mutate(X_00010_00003 = ifelse(is.na(X_00010_00003), X_.Discontinued.Dec..31..2011._00010_00003, X_00010_00003)) %>%
  mutate(date = as.Date(Date),
         month = month(Date),
         day = day(Date),
         year = year(Date)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(X_00010_00003, na.rm = TRUE)) 
#should return equal nrows as first pull, since it's already daily
write_csv(temp.dv, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_", gsub("-", "", Sys.Date()), ".csv"))

#note! uv is larger mem; almost 10 GB ram with other stuff loaded
temp.uv <- readNWISuv(siteNumbers = unique(g.get[g.get$data_type_cd == "uv", ]$site_no), 
                      parameterCd = "00010",
                      startDate = min(g.get[g.get$data_type_cd == "uv", ]$begin_date),
                      endDate = max(g.get[g.get$data_type_cd == "uv", ]$end_date))
temp.uv <- temp.uv %>%
  mutate(X_00010_00000 = ifelse(is.na(X_00010_00000), X_Discontinued.Sept..24..2014...Discontinued.Sept..24..2014._00010_00000, X_00010_00000)) %>%
  mutate(X_00010_00000 = ifelse(is.na(X_00010_00000), X_.YSI.EXO._00010_00000, X_00010_00000)) %>%
  mutate(date = as.Date(dateTime),
         month = month(dateTime),
         day = day(dateTime),
         year = year(dateTime)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(X_00010_00000, na.rm = T)) %>% ungroup()
write_csv(temp.uv, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_uvonly_", gsub("-", "", Sys.Date()), ".csv"))

temp.qw <- readNWISqw(siteNumbers = unique(g.get[g.get$data_type_cd == "qw", ]$site_no), 
                      parameterCd = "00010",
                      startDate = min(g.get[g.get$data_type_cd == "qw", ]$begin_date),
                      endDate = max(g.get[g.get$data_type_cd == "qw", ]$end_date))
temp.qw <- temp.qw[is.na(temp.qw$sample_end_dt),] #get rid of representing a sample period instead of a day (seems to be more than one day?)
temp.qw <- temp.qw %>%
  mutate(date = as.Date(sample_dt),
         month = month(date),
         day = day(date),
         year = year(date)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(result_va, na.rm = T)) %>% ungroup()
write_csv(temp.qw, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_qwonly_", gsub("-", "", Sys.Date()), ".csv"))

#combine all three
temp.dv <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_dvonly_20231024.csv")) %>% filter(!is.na(mn_daily_temp))
temp.qw <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_qwonly_20231024.csv")) %>% filter(!is.na(mn_daily_temp))
temp.uv <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_ARMOOK_uvonly_20231025.csv")) %>% filter(!is.na(mn_daily_temp))
#remove repeat data by level
tmp <- temp.uv %>% anti_join(temp.dv %>% select(site_no, year, month, day))
temp.all <- bind_rows(temp.dv, tmp)
tmp <- temp.qw %>% anti_join(temp.all %>% select(site_no, year, month, day)) #get qw data not already in combined data + add in
temp.all <- bind_rows(temp.all, tmp)

write_csv(temp.all, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_daily_temp_ARMOOK_all_", gsub("-", "", Sys.Date()), ".csv"))

#2. Calculate monthly mean temperature for each month of data ----
wtemp <- temp.all %>%
  group_by(site_no, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp, na.rm = T),
            ct = n()) %>%
  arrange(year, month, site_no)

write_csv(wtemp, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_monthly_temp_ARMOOK_all_", gsub("-", "", Sys.Date()), ".csv"))

rm(temp.dv, temp.qw, temp.uv, temp.all, tmp)

#site locations
xy <- g.info %>% 
  select(site_no, station_nm, dec_lat_va, dec_long_va) %>% 
  filter(site_no %in% wtemp$site_no) %>% 
  distinct()
write_csv(xy, paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_temp_locations.csv"))

#2+. Clean, compile, collected temperature data ----
dirs <- list.dirs(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data"), recursive = FALSE)

##eel ----
h <- readxl::read_excel(list.files(path = dirs[[1]], pattern = "[Hh]eader", full.names = TRUE), skip = 6) #get site info
t <- lapply(list.files(dirs[[1]], "[0-9].csv", full.names = TRUE), read_csv, skip = 1) #read in individual temp files
#one was wonky, so I corrected the column names in command line and reloaded 
t[[8]] <- read_delim(list.files(dirs[[1]], "[0-9].csv", full.names = TRUE)[[8]], delim = "\t", escape_double = FALSE) 

for(i in 1:length(t)){
  colnames(t[[i]]) <- c("row_id", "date", "temp") #rename column for consistency
  t[[i]]$FileName <- gsub(".csv", "", list.files(dirs[[1]], ".csv")[[i]])
  
  t[[i]] <- t[[i]] %>% 
    mutate(date = sub(" .*$", "", date),
           n.c = nchar(date))
  if(all(unique(t[[i]]$n.c) == 8)) { #some files have diff date format
    t[[i]] <- t[[i]] %>%
      mutate(date = as.Date(date, format = "%m/%d/%y"),
             year = year(date),
             month = month(date),
             day = day(date)) %>%
      select(-n.c)
  } else {
    t[[i]] <- t[[i]] %>%
      mutate(date = as.Date(date, format = "%m/%d/%Y"),
             year = year(date),
             month = month(date),
             day = day(date)) %>%
      select(-n.c)
  }
}

t <- bind_rows(t)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Eel_Header__Data/eel_water_temp_all_comped.csv"))

t <- t %>% 
  select(-row_id) %>%
  group_by(FileName, year, month, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% ungroup()

t <- t %>%
  group_by(FileName, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp, na.rm = T),
            ct = n()) %>%
  arrange(year, month, FileName)
t <- t %>%
  left_join(., h %>% select(FileName, Lat, Long))
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/eel_water_monthly_temp.csv"))

##MO JODI----
h <- read_csv(list.files(path = dirs[[2]], pattern = "[Hh]eader", full.names = TRUE))
t <- lapply(list.files(dirs[[2]], "[0-9].csv", full.names = TRUE), read_csv) #read in individual temp files

t <- lapply(t, function(x) {
  df <- x %>% select(SiteID, Date2, Temp_C)
  colnames(df) <- c("site_no", "date", "temp")
  df <- df %>% 
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           year = year(date),
           month = month(date),
           day = day(date)) %>%
    group_by(site_no, year, month, day) %>%
    summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% ungroup()
  return(df)
})

t <- bind_rows(t)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/MO_Jodi/MO_water_daily_temp.csv"))

t <- t %>%
  group_by(site_no, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp, na.rm = T),
            ct = n()) %>%
  arrange(year, month, site_no) %>% ungroup()
t <- t %>%
  left_join(., h %>% select(SiteID, Lat, Long), by = c("site_no" = "SiteID"))
#NOTE, SOME SITES DON'T HAVE LAT LONG DATA ??? FOR SOME REASON!
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/MO_water_monthly_temp.csv"))

## strawberry ----
h <- readxl::read_excel(list.files(path = dirs[[3]], pattern = "[Hh]eader", full.names = TRUE), skip = 6)
t <- lapply(list.files(dirs[[3]], "[0-9].xlsx", full.names = TRUE), readxl::read_excel) #read in individual temp files

for(i in 1:length(t)){
  colnames(t[[i]]) <- c("date", "temp") #rename column for consistency
  t[[i]]$FileName <- gsub(".xlsx", "", list.files(dirs[[3]], "[0-9].xlsx")[[i]])
  
  t[[i]] <- t[[i]] %>%
    mutate(date = sub(" .*$", "", date),
           date = as.Date(date, format = "%Y-%m-%d"),
           year = year(date),
           month = month(date),
           day = day(date))
}

t <- bind_rows(t)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/Strawberry_Header__Data/strawb_water_temp_all_comped.csv"))

#daily temp
t <- t %>% 
  group_by(FileName, year, month, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% ungroup()

#monthly temp
t <- t %>%
  group_by(FileName, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp, na.rm = T),
            ct = n()) %>%
  arrange(year, month, FileName)
t <- t %>%
  left_join(., h %>% select(FileName, Lat, Long)) #add in lat long
#missed a misspelling
fix <- t %>%
  filter(is.na(Lat)) %>%
  mutate(FileName = ifelse(FileName == "10996600__Straw_HulettR_9.23.2019", "10996600__Straw_HulettRd_9.23.2019", FileName))
fix <- fix %>% select(-Lat, -Long) %>% left_join(., h %>% select(FileName, Lat, Long))

t <- bind_rows(t %>% filter(!is.na(Lat)), fix)

write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/strawb_water_monthly_temp.csv"))

##adams lab ----
h <- read_csv(list.files(path = dirs[[4]], pattern = "[Hh]eader", full.names = TRUE), skip = 6)
t <- lapply(list.files(dirs[[4]], ".csv", full.names = TRUE)[-22], read_csv) #read in individual temp files, ignore header file

for(i in 1:length(t)){
  t[[i]] <- t[[i]] %>% select(Date, Temperature_C) %>% rename(date = Date, temp = Temperature_C) #rename column for consistency
  t[[i]]$FileName <- gsub(".csv", "", list.files(dirs[[4]], ".csv")[-22][[i]])
  
  t[[i]] <- t[[i]] %>%
    mutate(date = as.Date(date, format = "%Y-%m-%d"),
           year = year(date),
           month = month(date),
           day = day(date))
}

t <- bind_rows(t)
write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/UCA_AdamsLab/uca_water_temp_all_comped.csv"))

#daily temp
t <- t %>% 
  group_by(FileName, year, month, day) %>%
  summarise(mn_daily_temp = mean(temp, na.rm = T)) %>% ungroup()

#monthly temp
t <- t %>%
  group_by(FileName, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp, na.rm = T),
            ct = n()) %>%
  arrange(year, month, FileName)
t <- t %>%
  left_join(., h %>% select(FileName, Lat, Long) %>% mutate(FileName = gsub(".csv", "", FileName))) #add in lat long

write_csv(t, paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data/uca_water_monthly_temp.csv"))
rm(h, t)

##compile ALL data into single file, with location info ----
nu <- lapply(list.files(paste0(PATH, "/02_EnvDat/raw_stream_temp/non_usgs_temp_data"), ".csv", full.names = TRUE), read_csv)
us <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/USGS_water_monthly_temp_ARMOOK_all_20231025.csv"))
us <- left_join(us, xy %>% rename(Lat = dec_lat_va, Long = dec_long_va) %>% select(-station_nm))

nu <- lapply(nu, function(x) { 
  if(any(grepl("FileName", colnames(x)))) {
    x <- x %>% rename(site_no = FileName)
  }
  
  x <- x %>% mutate(site_no = as.character(site_no))
  
  return(x)
})

nu <- bind_rows(nu)

wtemp <- bind_rows(nu, us)
write_csv(wtemp, paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_temp_", gsub("-", "", Sys.Date()), ".csv"))

rm(nu, us)

#3. Exclude months with fewer than 20 days temp data ----
wtemp <- wtemp %>% ungroup() %>% filter(ct >= 20) %>%
  mutate(date = with(., sprintf("%d-%02d", year, month))) %>%
  filter(!is.na(Lat))

length(unique(wtemp$site_no)) #only 49 of the original USGS gauge stations have mean monthly water temperature data with > 20 records, but 306 sites total

#4. Collect monthly mean air temp from the National Oceanic and Atmospheric Agency National Center for Climate Data ----
##4a. From the county where each site was located for a corresponding time period to the water temp data
#NOAA county location IDs:
ids <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/county_fips_master.csv")) %>%
  mutate(fips = case_when(nchar(as.character(fips)) == 4 ~ paste0("0", as.character(fips)),
                          T ~ as.character(fips))) %>%
  filter(state_abbr %in% c("AR", "OK", "MO"))

#find which counties the temp locations are in:
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)

sf_use_s2(FALSE) #had to turn off spherical geometry to do this, shouldn't be an issue (??)

xy <- st_as_sf(wtemp %>% select(site_no, Long, Lat) %>% distinct(), coords = c("Long", "Lat"), crs = 4269)

#find which counties each site falls in
str.noaa <- hlnd %>% 
  st_join(., xy) %>% 
  filter(!is.na(site_no))

counties <- wtemp %>%
  left_join(., str.noaa %>% st_drop_geometry()) %>% #add back in county names, states
  group_by(ID) %>%
  reframe(start = min(date), end = max(date), months = n()) %>% #get search date range for each county with temp sites
  rename(county = ID) %>%
  distinct() %>%
  mutate(state_abbr = case_when(grepl("arkansas", county, fixed = TRUE) ~ "AR",
                                grepl("missouri", county, fixed = TRUE) ~ "MO",
                                grepl("oklahoma", county, fixed = TRUE) ~ "OK",
                                T ~ NA),
         county_name = sub("^.*?,", "", county),
         county_name = paste(str_to_title(county_name), "County"),
         county_name = case_when(grepl("Mc", county_name, fixed = TRUE) ~ gsub("(Mc)([a-z])", "\\1\\U\\2", county_name, perl = TRUE),
                                 T ~ county_name)) %>%
  left_join(., ids %>% select(county_name, state_abbr, fips), by = c("county_name", "state_abbr")) %>% #add in fips code for search
  #noaa search requires days, so get start and end of each month in search range
  mutate(p_start = as.Date(paste0(start, "-01")), 
         p_end = ceiling_date(ymd(paste0(end, "-01")), 'month') - days(1),
         n_days = difftime(p_end, p_start))

#get daily temperature data for each county
base_url <- "https://www.ncei.noaa.gov/cdo-web/api/v2/data"

##get dates in chunks for each county in 12 month increments (seems to be limit?), make a dataframe for search to pull from
d_df <- data.frame()
for(j in 1:nrow(counties)) {
  tmp <- counties[j, ]
  df12 <- data.frame(fips = tmp$fips, 
                     start = seq(tmp$p_start, tmp$p_end, by = "13 months"), 
                     end = NA, 
                     status = NA)
  df12 <- df12 %>% 
    mutate(end = lead(start) - days(1), 
           end = case_when(is.na(end) ~ tmp$p_end,
                           T ~ end))
  d_df <- bind_rows(d_df, df12)
}

a_temp <- data.frame()
for (i in 1:nrow(d_df)) {
  #make get url
  url <- paste0(base_url,
                "?datasetid=GHCND",
                "&datatypeid=TMAX&datatypeid=TMIN",
                "&locationid=FIPS:", d_df[i, "fips"],
                "&startdate=", d_df[i, "start"],
                "&enddate=", d_df[i, "end"],
                "&units=metric",
                "&limit=1000")
  req <- request(url) %>%
    req_headers(token = api_token) %>%
    req_error(is_error = function(resp) FALSE) %>% #don't stop if it's a bad request
    req_perform()
  d_df[i, "status"] <- req$status_code #to check specific sections that go wrong
  
  if(req$status_code == 200) {
    
    res <- resp_body_json(req, simplifyVector = TRUE)$results
    
    if(all(is.na(res$value))) { #for searches that go through but don't return any data
      
      d_df[i, "status"] <- NA
      cat(paste0(d_df[i, "fips"], ": ", d_df[i, "start"], " to ", d_df[i, "end"]), "all NA..\n")
      
    } else {
      #flag returns greater than given limit
      if(resp_body_json(req, simplifyVector = TRUE)$metadata$resultset$count > 1000) { 
        d_df[i, "status"] <- resp_body_json(req, simplifyVector = TRUE)$metadata$resultset$count
      }
      
      res$fips <- d_df[i, "fips"]
      a_temp <- bind_rows(a_temp, res)
      cat(paste0(d_df[i, "fips"], ": ", d_df[i, "start"], " to ", d_df[i, "end"]), "complete..\n")
      
    }
    
  } else {
    
    cat(paste0(d_df[i, "fips"], ": ",d_df[i, "start"], " to ", d_df[i, "end"]), "bad request..\n")
    
  }
  
  if(i == nrow(d_df)) { cat("\nDone!\n") }
}

d_df %>% filter(status != 200 | is.na(status)) #check for bad or NA requests and redo
write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv")) #save just in case

#redo the bad requests...
d_df2 <- d_df %>% filter(status == 503)
a_temp2 <- data.frame()
for (i in 1:nrow(d_df2)) {
  #make get url
  url <- paste0(base_url,
                "?datasetid=GHCND",
                "&datatypeid=TMAX&datatypeid=TMIN",
                "&locationid=FIPS:", d_df2[i, "fips"],
                "&startdate=", d_df2[i, "start"],
                "&enddate=", d_df2[i, "end"],
                "&units=metric",
                "&limit=1000")
  req <- request(url) %>%
    req_headers(token = api_token) %>%
    req_error(is_error = function(resp) FALSE) %>% #don't stop if it's a bad request
    req_perform()
  d_df2[i, "status"] <- req$status_code #to check specific sections that go wrong
  
  if(req$status_code == 200) {
    
    res <- resp_body_json(req, simplifyVector = TRUE)$results
    
    if(all(is.na(res$value))) { #for searches that go through but don't return any data
      
      d_df2[i, "status"] <- NA
      cat(paste0(d_df2[i, "fips"], ": ", d_df2[i, "start"], " to ", d_df2[i, "end"]), "all NA..\n")
      
    } else {
      #flag returns greater than given limit
      if(resp_body_json(req, simplifyVector = TRUE)$metadata$resultset$count > 1000) { 
        d_df2[i, "status"] <- resp_body_json(req, simplifyVector = TRUE)$metadata$resultset$count
      }
      
      res$fips <- d_df2[i, "fips"]
      a_temp2 <- bind_rows(a_temp2, res)
      cat(paste0(d_df2[i, "fips"], ": ", d_df2[i, "start"], " to ", d_df2[i, "end"]), "complete..\n")
      
    }
    
  } else {
    
    cat(paste0(d_df2[i, "fips"], ": ",d_df2[i, "start"], " to ", d_df2[i, "end"]), "bad request..\n")
    
  }
  
  if(i == nrow(d_df2)) { cat("\nDone!\n") }
}

a_temp <- bind_rows(a_temp, a_temp2)
write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv")) #save just in case

d_df <- bind_rows(d_df %>% anti_join(d_df2, by = c("fips", "start", "end")), d_df2) #update statuses, then rerun d_df2 search at L409

#check the NA runs... check: https://www.ncdc.noaa.gov/cdo-web/search , it's likely that temp is not available at those queries
rm(req, tmp, df12, res)

d_df2 <- d_df %>% filter(status != 200)
d_df2$status_2 <- NA
a_temp2 <- data.frame()
for(i in 1:nrow(d_df2)) {
  n <- d_df2[i, "status"]
  
  tot <- seq(1000, n, by = 1000)
  
  sub_res <- data.frame()
  for(s in 1:length(tot)) {
    url <- paste0(base_url,
                  "?datasetid=GHCND",
                  "&datatypeid=TMAX&datatypeid=TMIN",
                  "&locationid=FIPS:", d_df2[i, "fips"],
                  "&startdate=", d_df2[i, "start"],
                  "&enddate=", d_df2[i, "end"],
                  "&units=metric",
                  "&limit=1000", 
                  "&offset=", tot[s])
    req <- request(url) %>%
      req_headers(token = api_token) %>%
      req_error(is_error = function(resp) FALSE) %>% #don't stop if it's a bad request
      req_perform()
    
    if(req$status_code == 200) {
      res <- resp_body_json(req, simplifyVector = TRUE)$results
      sub_res <- bind_rows(sub_res, res)
      cat(paste0(d_df2[i, "fips"], ": ", d_df2[i, "start"], " to ", d_df2[i, "end"]), s, "offset complete..\n")
    } else {
      d_df2[i, "status_2"] <- "redo"
      cat(paste0(d_df2[i, "fips"], ": ", d_df2[i, "start"], " to ", d_df2[i, "end"]), s, "offset bad request..\n")
    }
    
    if(s == length(tot)) { cat(paste0(d_df2[i, "fips"], ": ", d_df2[i, "start"], " to ", d_df2[i, "end"]), "done!\n") }
  }
  
  if(nrow(sub_res) != 0) {
    sub_res$fips <- d_df2[i, "fips"]
    a_temp2 <- bind_rows(a_temp2, sub_res)
  }
}

#redo request marked redo:
d_df3 <- d_df2 %>% filter(status_2 == "redo")
a_temp <- bind_rows(a_temp, a_temp2) #add back in the data that did return and save
write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv")) #save just in case

a_temp2 <- data.frame()
for(i in 1:nrow(d_df3)) {
  n <- d_df3[i, "status"]
  
  tot <- seq(1000, n, by = 1000)
  
  sub_res <- data.frame()
  for(s in 1:length(tot)) {
    url <- paste0(base_url,
                  "?datasetid=GHCND",
                  "&datatypeid=TMAX&datatypeid=TMIN",
                  "&locationid=FIPS:", d_df3[i, "fips"],
                  "&startdate=", d_df3[i, "start"],
                  "&enddate=", d_df3[i, "end"],
                  "&units=metric",
                  "&limit=1000", 
                  "&offset=", tot[s])
    req <- request(url) %>%
      req_headers(token = api_token) %>%
      req_error(is_error = function(resp) FALSE) %>% #don't stop if it's a bad request
      req_perform()
    
    if(req$status_code == 200) {
      d_df3[i, "status_2"] <- NA
      res <- resp_body_json(req, simplifyVector = TRUE)$results
      sub_res <- bind_rows(sub_res, res)
      cat(paste0(d_df3[i, "fips"], ": ", d_df3[i, "start"], " to ", d_df3[i, "end"]), s, "offset complete..\n")
    } else {
      d_df3[i, "status_2"] <- "redo"
      cat(paste0(d_df3[i, "fips"], ": ", d_df3[i, "start"], " to ", d_df3[i, "end"]), s, "offset bad request..\n")
    }
    
    if(s == length(tot)) { cat(paste0(d_df3[i, "fips"], ": ", d_df3[i, "start"], " to ", d_df3[i, "end"]), "done!\n") }
  }
  
  if(nrow(sub_res) != 0) {
    sub_res$fips <- d_df3[i, "fips"]
    a_temp2 <- bind_rows(a_temp2, sub_res)
  }
}

d_df2 <- bind_rows(d_df2 %>% anti_join(d_df3, by = c("fips", "start", "end")), d_df3) #update statuses, then rerun d_df3 search at L503
a_temp <- bind_rows(a_temp, a_temp2) #add back in the data that did return and save
write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv")) #save just in case

d_df3 <- d_df2 %>% filter(status_2 == "redo")

a_temp <- a_temp %>% distinct() #get rid of any accidental duplicates
write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv")) #save just in case

rm(a_temp2, req, sub_res, res, i, j, n, s, tot, d_df2, d_df3)

#calc monthly means
a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_daily_temp_countylevel_20231025.csv"))
a_mn_temp <- a_temp %>%
  mutate(year = year(date),
         month = month(date)) %>%
  select(-attributes) %>% 
  pivot_wider(names_from = datatype, values_from = value) %>%
  group_by(fips, year, month) %>%
  summarise(mn_max_monthly_temp = mean(TMAX, na.rm = T),
            mn_min_monthly_temp = mean(TMIN, na.rm = T)) %>%
  rowwise() %>%
  mutate(mn_monthly_temp = mean(c(mn_max_monthly_temp, mn_min_monthly_temp), na.rm = T)) %>%
  ungroup() %>%
  mutate(date = with(., sprintf("%d-%02d", year, month)),
         p_date = as.Date(paste0(date, "-01"))) %>%
  arrange(fips, p_date)
write_csv(a_mn_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_monthly_temp_countylevel_", gsub("-", "", Sys.Date()), ".csv"))

a_mn_temp_id <- a_mn_temp %>%
  left_join(., counties %>% select(fips, county)) %>%
  left_join(., str.noaa %>% st_drop_geometry(), 
            by = c("county" = "ID"), multiple = "all")
write_csv(a_mn_temp_id, paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_monthly_temp_countylevel_waterIDS.csv"))

rm(a_mn_temp, a_mn_temp_id, a_temp, d_df, ids, str.noaa, counties, wtemp, xy)

#5. Attach gage data to streams + add flow cat ----
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% 
  st_transform(5070)
temp.sites <- st_as_sf(w_temp %>% select(site_no, Lat, Long) %>% filter(!is.na(Lat)) %>% distinct(), coords = c("Long", "Lat"), crs = 4269)
alb.sites <- temp.sites %>% st_transform(5070)

#assign COMID to each site + calc distance to nearest line (meters)
nr.line <- st_nearest_feature(alb.sites, nhd, check_crs = TRUE) #find nearest nhd strm
nhd_near <- nhd[nr.line,]
temp.sites$COMID <- nhd_near$COMID #add comid match to sites
dist <- as.vector(st_distance(alb.sites, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm
temp.sites$dist2strm_m_nhd <- dist

#assign flow_type to each site + calc distance to nearest line (meters)
flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(5070)
flw <- flw %>% filter(Flow_type != "BigR")
nr.line <- st_nearest_feature(alb.sites, flw, check_crs = TRUE) #find nearest strm not river
flw_near <- flw[nr.line,]
temp.sites$flw_name <- flw_near$Name #add strm names
temp.sites$flw_type <- flw_near$Flow_type #add flow types
dist <- as.vector(st_distance(alb.sites, flw_near, by_element = TRUE)) #get distance from site to strm assigned
temp.sites$dist2strm_m_flw <- dist

rm(nhd, flw, alb.sites, flw_near, nhd_near, dist, nr.line)

#save location info
temp.sites <- temp.sites %>%
  mutate(Long = st_coordinates(.)[,"X"],
         Lat = st_coordinates(.)[,"Y"],
         across(contains("dist2strm"), ~ round(.x, digits = 3)))
write_csv(temp.sites, paste0(PATH, "/02_EnvDat/raw_stream_temp/all_water_temp_locations_info.csv"))

#6. Get USGS climate data ----
w_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_temp_20231025.csv")) %>%
  rename(water_temp = mn_monthly_temp) %>%
  filter(ct >= 20) %>%
  mutate(date = with(., sprintf("%d-%02d", year, month)),
         p_date = as.Date(paste0(date, "-01")))

#Find which counties each site falls in again, also all gage sites for hit dataset
ids <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/county_fips_master.csv")) %>%
  mutate(fips = case_when(nchar(as.character(fips)) == 4 ~ paste0("0", as.character(fips)),
                          T ~ as.character(fips))) %>%
  filter(state_abbr %in% c("AR", "OK", "MO"))

#get site numbers from data:
g.nums <- unique(hit$STAID_0)

#get location info for gages
g.info <- whatNWISdata(siteNumbers = g.nums) %>% 
  select(site_no, dec_lat_va, dec_long_va) %>%
  distinct() %>%
  rename(Lat = dec_lat_va,
         Long = dec_long_va)
other.dat <- w_temp %>% select(site_no, Lat, Long) %>% filter(!is.na(Lat)) %>% distinct() %>% anti_join(g.info)

g.xy <- bind_rows(g.info, other.dat) %>%
  st_as_sf(., coords = c("Long", "Lat"), crs = 4269) #combine usgs + other data

#get counties in ARMOOK
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)

sf_use_s2(FALSE)

counties <- hlnd %>% 
  st_join(., g.xy) %>% 
  filter(!is.na(site_no)) %>% 
  rename(county = ID) %>%
  distinct() %>%
  mutate(state_abbr = case_when(grepl("arkansas", county, fixed = TRUE) ~ "AR",
                                grepl("missouri", county, fixed = TRUE) ~ "MO",
                                grepl("oklahoma", county, fixed = TRUE) ~ "OK",
                                T ~ NA),
         county_name = sub("^.*?,", "", county),
         county_name = paste(str_to_title(county_name), "County"),
         county_name = ifelse(grepl("Mc", county_name, fixed = TRUE), gsub("(Mc)([a-z])", "\\1\\U\\2", county_name, perl = TRUE), county_name),
         county_name = ifelse(grepl("St ", county_name, fixed = TRUE), gsub("St ", "St. ", county_name), county_name), 
         county_name = ifelse(grepl("St. Louis City", county_name, fixed = TRUE), "St. Louis city", county_name)) %>%
  left_join(., ids %>% select(county_name, state_abbr, fips), by = c("county_name", "state_abbr")) #add in fips code for search
# write_csv(counties, paste0(PATH, "/02_EnvDat/usgs_gage_prevHITsites_countyinfo.csv"))

#download by county
uniq.fips <- unique(counties$fips)
for(i in 1:2) {
  fips <- uniq.fips[i]
  cat("Starting", fips, "county download...\n")
  # fips <- "05019"
  save.dest <- paste0(PATH, "/02_EnvDat/raw_air_temp/usgs_nccv_rawdata/USGS_climate_", fips, ".csv")
  
  if(file.exists(save.dest)) {
    cat(fips, "county data already exists.\n"); next
  } else {
    #set download URL
    url <- paste0("https://regclim.ceoas.oregonstate.edu/nccv2/maca2/county/csv//US",
                  fips, "/US", fips, "_MeanModel_metric.csv")
    user_agent <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36"
    
    out <- suppressWarnings(
      try(
        download.file(url = url, 
                      destfile = save.dest, 
                      quiet = TRUE, 
                      headers = c("User-Agent" = user_agent)),
        silent = TRUE))
    
    if(inherits(out, "try-error")) {
      cat(fips, "URL not found..\n")
    } else {
      cat(fips, "county downloaded..\n")
    }
    
    Sys.sleep(2) #wait!!!
    
  }
}

#clean usgs data for use 
t.clim <- lapply(list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/usgs_nccv_rawdata/"), pattern = ".csv", full.names = TRUE), read_csv, skip = 20)
fips <- sub(".csv", "", sub("USGS_climate_", "", list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/usgs_nccv_rawdata/"), pattern = ".csv")))

for(i in 1:length(t.clim)) {
  t.clim[[i]] <- t.clim[[i]] %>%
    select(Date, contains("temp")) %>%
    rename(date = Date,
           mn_temp_4.5 = `RCP4.5 Mean temperature (deg_C)`,
           max_temp_4.5 = `RCP4.5 Max temperature (deg_C)`,
           min_temp_4.5 = `RCP4.5 Min temperature (deg_C)`,
           mn_temp_8.5 = `RCP8.5 Mean temperature (deg_C)`,
           max_temp_8.5 = `RCP8.5 Max temperature (deg_C)`,
           min_temp_8.5 = `RCP8.5 Min temperature (deg_C)`) %>%
    mutate(fips_code = fips[[i]],
           date = as.Date(date, format = "%m/%d/%Y"))
}

t.clim <- bind_rows(t.clim)

write_csv(t.clim, paste0(PATH, "/02_EnvDat/raw_air_temp/usgs_air_monthly_temp_countylevel.csv"))

# + Plot data available ----
#armook counties
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)
#county fips
ids <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/county_fips_master.csv")) %>%
  mutate(fips = case_when(nchar(as.character(fips)) == 4 ~ paste0("0", as.character(fips)),
                          T ~ as.character(fips))) %>%
  filter(state_abbr %in% c("AR", "OK", "MO"))
#ecoregion with flow types
eco <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/ecoreg_l3_interior_highlands_crop.shp")) %>%
  st_union()
water_sites <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/all_water_temp_locations_info.csv")) %>%
  mutate(source = case_when(grepl("^[0-9]+$", site_no) & nchar(site_no) >= 8 ~ "USGS Gage",
                            grepl("^[0-9]+$", site_no) & nchar(site_no) <= 3 ~ "MO_Jodi",
                            grepl("^1099", site_no) & grepl("[Ss]traw", site_no) ~ "Strawberry",
                            grepl("^[0-9]+", site_no) & grepl("Ouachita", site_no) ~ "Eel",
                            T ~ "UCA_AdamsLab")) %>%
  st_as_sf(., coords = c("Long", "Lat"), crs = 4269)
#county shp to fips join
# sf_use_s2(FALSE)
eco.buff <- eco %>% st_buffer(dist = 0.09)

water_sites <- water_sites %>% 
  mutate(in_eco_reg = st_within(., eco.buff),
         in_eco_reg = ifelse(as.numeric(in_eco_reg) == 1, TRUE, FALSE),
         in_eco_reg = ifelse(is.na(in_eco_reg), FALSE, in_eco_reg))

counties <- hlnd %>% 
  st_join(., water_sites) %>% 
  filter(!is.na(site_no)) %>% 
  rename(county = ID) %>%
  distinct() %>%
  mutate(state_abbr = case_when(grepl("arkansas", county, fixed = TRUE) ~ "AR",
                                grepl("missouri", county, fixed = TRUE) ~ "MO",
                                grepl("oklahoma", county, fixed = TRUE) ~ "OK",
                                T ~ NA),
         county_name = sub("^.*?,", "", county),
         county_name = paste(str_to_title(county_name), "County"),
         county_name = ifelse(grepl("Mc", county_name, fixed = TRUE), gsub("(Mc)([a-z])", "\\1\\U\\2", county_name, perl = TRUE), county_name),
         county_name = ifelse(grepl("St ", county_name, fixed = TRUE), gsub("St ", "St. ", county_name), county_name), 
         county_name = ifelse(grepl("St. Louis City", county_name, fixed = TRUE), "St. Louis city", county_name)) %>%
  left_join(., ids %>% select(county_name, state_abbr, fips), by = c("county_name", "state_abbr")) #add in fips code for search
# write_csv(counties %>% st_drop_geometry(), paste0(PATH, "/02_EnvDat/raw_stream_temp/all_water_temp_locations_countyflowinfo.csv"))

p1 <- ggplot() +
  geom_sf(data = hlnd) +
  geom_sf(data = counties, fill = "darkgrey") +
  geom_sf(data = eco, fill = NA, color = "black", lwd = 1) +
  geom_sf(data = water_sites, aes(shape = source, color = flw_type), size = 2) +
  scale_color_viridis_d(end = 0.9) +
  theme_minimal()
# p1
ggsave(paste0(PATH, "/99_figures/fine_scale_temp_site_map.png"), plot = p1, dpi = 300, width = 8, height = 6, bg = "white")





