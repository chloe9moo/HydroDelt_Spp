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
library(dataRetrieval)

PATH <- getwd()

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
#find which counties the USGS stations are in:
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)
counties <- hlnd %>% 
  st_join(., st_as_sf(xy, coords = c("dec_long_va", "dec_lat_va"), crs = 4269)) %>% #had to turn off spherical geometry to do this
  filter(!is.na(site_no))
temp %>% 
  left_join(., counties %>% st_drop_geometry()) %>%
  group_by(site_no) %>%
  reframe(start = min(date), end = max(date), county = ID, months = n()) %>%
  distinct()
#NOAA location IDs:
ids <- read_csv(paste0(PATH, "/01_Data/county_fips_master.csv")) %>% 
  mutate(fips = paste0("0", as.character(fips))) %>%
  filter(state_abbr %in% c("AR", "OK", "MO"))


clim <- read_csv(list.files(paste0(PATH, "/01_Data/daily_climate/"), full.names = T)) %>%
  filter(!is.na(TOBS)) %>%
  mutate(year = year(DATE),
         month = month(DATE)) %>%
  group_by(year, month) %>%
  summarise(mn_monthly_temp = mean(TOBS)) %>%
  ungroup() %>%
  mutate(date = as.Date(with(., sprintf("%d-%02d", year, month)))) %>%
  arrange(date)
#this data is from izard county as a test! also date is not recognized (can't plot) likely because it doesn't have a day
#also, I think I need min temp and max temp, and then take a mean of them all (see paper)




ggplot() +
  geom_sf(data = hlnd) +
  geom_sf(data = counties, fill = "red") +
  geom_sf(data = st_as_sf(xy, coords = c("dec_long_va", "dec_lat_va"), crs = 4269)) +
  theme_minimal()


