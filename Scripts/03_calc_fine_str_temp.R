## Prep + calculate fine scale stream temperature
# by C.E.Moore V.
# 
# following the methods of Middaugh et al. 2016 (Ecology of Freshwater Fish)

## STEPS:
#1. Download temperature data for USGS gages
#2. Calculate monthly mean temperature for each month of data
#3. Exclude months with fewer than 20 days temp data
#4. Collect monthly mean air temp from the National Oceanic and Atmospheric Agency National Center for Climate Data
##4a. From the county where each site was located for a corresponding time period to the water temp data
#5. Calculate least-squares linear regression for each site to predict stream temperature from air temperature
#6. Using the USGS National Climate Change Viewer to get historical and future air temperatures (not sure if this is necessary?)

library(sf); library(tidyverse); library(maps)
library(dataRetrieval); library(httr2)

PATH <- getwd()

source(paste0(PATH, "/Scripts/XX_api_token.R"))

fish.dat <- read_csv(list.files(paste0(PATH, "/01_Data"), "Fishes", full.name=TRUE))
bug.dat <- read_csv(list.files(paste0(PATH, "/01_Data"), "Insect", full.name=TRUE))

#1. Download temperature data for USGS gages ----
# https://waterdata.usgs.gov/blog/dataretrieval/
#get site numbers from data:
g.nums <- unique(c(unique(fish.dat$STAID), unique(bug.dat$STAID)))
g.nums <- as.character(g.nums)
g.nums <- unlist(lapply(g.nums, function(x) ifelse(nchar(x) < 8, paste0("0", x), x)))

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
  mutate(date = as.Date(Date),
         month = month(Date),
         day = day(Date),
         year = year(Date)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(X_00010_00003))

#note! uv is larger mem; almost 5 GB ram itself
temp.uv <- readNWISuv(siteNumbers = unique(g.get[g.get$data_type_cd == "uv", ]$site_no), 
                      parameterCd = "00010",
                      startDate = min(g.get[g.get$data_type_cd == "uv", ]$begin_date),
                      endDate = max(g.get[g.get$data_type_cd == "uv", ]$end_date))
temp.uv <- temp.uv %>%
  mutate(date = as.Date(dateTime),
         month = month(dateTime),
         day = day(dateTime),
         year = year(dateTime)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(X_00010_00000)) %>% ungroup()

temp.qw <- readNWISqw(siteNumbers = unique(g.get[g.get$data_type_cd == "qw", ]$site_no), 
                      parameterCd = "00010",
                      startDate = min(g.get[g.get$data_type_cd == "qw", ]$begin_date),
                      endDate = max(g.get[g.get$data_type_cd == "qw", ]$end_date))
temp.qw <- temp.qw[is.na(temp.qw$sample_end_dt),] #get rid of apparent sample period rows (seems to be more than one day?)
temp.qw <- temp.qw %>%
  mutate(date = as.Date(sample_dt),
         month = month(sample_dt),
         day = day(sample_dt),
         year = year(sample_dt)) %>%
  group_by(site_no, year, month, day) %>%
  summarise(mn_daily_temp = mean(result_va)) %>% ungroup()

#2. Calculate monthly mean temperature for each month of data ----
temp <- bind_rows(temp.dv, temp.qw, temp.uv) %>%
  distinct(site_no, year, month, day, .keep_all = TRUE) %>%
  group_by(site_no, year, month) %>%
  summarise(mn_monthly_temp = mean(mn_daily_temp),
            ct = n()) %>%
  arrange(year, month, site_no)

#save just in case
bind_rows(temp.dv, temp.qw, temp.uv) %>%
  distinct(site_no, year, month, day, .keep_all = TRUE) %>%
  write_csv(., paste0(PATH, "/01_Data/USGS_gauge_daily_wtemp.csv"))
write_csv(temp, paste0(PATH, "/01_Data/USGS_gauge_monthly_wtemp.csv"))

rm(temp.dv, temp.qw, temp.uv)

#3. Exclude months with fewer than 20 days temp data ----
temp <- temp %>% filter(ct >= 20) %>%
  mutate(date = with(., sprintf("%d-%02d", year, month)))

length(unique(temp$site_no)) #only 23 of the original gauge stations have mean monthly water temperature data with > 20 records

#4. Collect monthly mean air temp from the National Oceanic and Atmospheric Agency National Center for Climate Data ----
##4a. From the county where each site was located for a corresponding time period to the water temp data

#site locations
xy <- g.info %>% 
  select(site_no, station_nm, dec_lat_va, dec_long_va) %>% 
  filter(site_no %in% temp$site_no) %>% 
  distinct()

#NOAA location IDs:
ids <- read_csv(paste0(PATH, "/01_Data/county_fips_master.csv")) %>%
  mutate(fips = case_when(nchar(as.character(fips)) == 4 ~ paste0("0", as.character(fips)),
                          T ~ as.character(fips))) %>%
  filter(state_abbr %in% c("AR", "OK", "MO"))

#find which counties the USGS stations are in:
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)
sf_use_s2(FALSE) #had to turn off spherical geometry to do this, shouldn't be an issue (??)
usgs.noaa <- hlnd %>% 
  st_join(., st_as_sf(xy, coords = c("dec_long_va", "dec_lat_va"), crs = 4269)) %>% 
  filter(!is.na(site_no))
counties <- temp %>%
  left_join(., usgs.noaa %>% st_drop_geometry()) %>%
  group_by(ID) %>%
  reframe(start = min(date), end = max(date), months = n()) %>%
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
  left_join(., ids %>% select(county_name, state_abbr, fips), by = c("county_name", "state_abbr")) %>%
  mutate(p_start = as.Date(paste0(start, "-01")),
         p_end = ceiling_date(ymd(paste0(end, "-01")), 'month') - days(1),
         n_days = difftime(p_end, p_start))

#get daily temperature data for each county
base_url <- "https://www.ncei.noaa.gov/cdo-web/api/v2/data"
##get dates in chunks for each county in 12 month increments (seems to be limit?)
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
# d_df <- d_df %>% filter(fips %in% redo$fips)
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
    res$fips <- d_df[i, "fips"]
    a_temp <- bind_rows(a_temp, res)
    cat(paste0(d_df[i, "fips"], ": ", d_df[i, "start"], " to ", d_df[i, "end"]), "complete..\n")
  } else {
    cat(paste0(d_df[i, "fips"], ": ",d_df[i, "start"], " to ", d_df[i, "end"]), "bad request..\n")
  }
  
  if(i == nrow(d_df)) { cat("\nDone!\n") }
}
rm(req, tmp, df12)

#check run
# a_temp %>% filter(is.na(value))

# write_csv(a_temp, paste0(PATH, "/01_Data/NOAA_daily_atemp.csv"))
# a_temp2 <- read_csv(paste0(PATH, "/01_Data/NOAA_daily_atemp.csv"))
# a_temp3 <- bind_rows(a_temp %>% mutate(date = as.Date(date)), a_temp2)
# a_temp <- a_temp3 %>% filter(!is.na(value))
write_csv(a_temp, paste0(PATH, "/01_Data/NOAA_daily_atemp.csv"))

#calc monthly means
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
write_csv(a_mn_temp, paste0(PATH, "/01_Data/NOAA_monthly_atemp.csv"))

a_mn_temp_id <- a_mn_temp %>%
  left_join(., counties %>% select(fips, county)) %>%
  left_join(., usgs.noaa %>% st_drop_geometry() %>% select(-station_nm), 
            by = c("county" = "ID"), multiple = "all")
write_csv(a_mn_temp_id, paste0(PATH, "/01_Data/NOAA_monthly_atemp_wIDS.csv"))

#5. Calculate least-squares linear regression for each site to predict stream temperature from air temperature ----
