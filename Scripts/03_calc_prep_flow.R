## FLOW CALC + ESTIMATION ##

library(tidyverse); library(sf)
library(dataRetrieval); library(EflowStats)
options(readr.show_col_types = FALSE)

PATH <- getwd()

# pull daily flow ----
min_yr_period <- 10 #set number of years req. for use in analyses

## get gage info + filter ----
flow.gage <- data.frame()
for(state_nm in c("AR", "MO", "OK")) {
  
  tmp <- whatNWISdata(stateCd = state_nm,
                      parameterCd = "00060", #discharge
                      service = "dv",
                      statCd = "00003")
  flow.gage <- rbind(flow.gage, tmp)
  
}
  
#save for later (JIC)
write_csv(flow.gage, paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"))

#filter gages for needs
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop") %>% #int. highlands boundary
  st_union() %>%
  st_buffer(., dist = 2000) #2 km buffer (for wiggle room)
# huc <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "huc10_interior_highlands_crop")

flow.gage <- flow.gage %>%
  mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
  #10 years of data
  filter(period >= min_yr_period * 365) %>%
  #make spatial for int high filter
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269) #NAD83
  
flow.gage <- flow.gage %>%
  mutate(int_high = as.logical(st_within(., eco)),
         int_high = ifelse(is.na(int_high), FALSE, int_high)) #use int. highland eco region + buffer

write_csv(st_drop_geometry(flow.gage), paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"))

#check
# ggplot() +
#   geom_sf(data = flow.gage, aes(color = int_high)) +
# # #   geom_sf(data = flow.gage, color = "blue") +
#   geom_sf(data = eco, fill = NA, color = "red") +
# # #   # geom_sf(data = huc, fill = NA, color = "blue") +
#   theme_minimal()

rm(tmp, eco, state_nm)

## retrieve flow data for filtered gages ----
flow.gage <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"), col_types = cols(site_no = col_character())) %>% #make sure sites load in correctly
  filter(int_high == TRUE)

site_nm <- unique(flow.gage$site_no)

for(i in seq_along(site_nm)) { #pull flow data and save
  nm <- site_nm[[i]]
  
  if(file.exists(paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))) { #in case it fails at any point
    message("Flow dat for ", nm, " already exists."); next
  } else {
    
    tmp.dat <- readNWISdv(siteNumbers = nm,
                          parameterCd = "00060",
                          statCd = "00003",
                          startDate = min(flow.gage[flow.gage$site_no == nm, ]$begin_date),
                          endDate = max(flow.gage[flow.gage$site_no == nm, ]$end_date))
    
    tmp.dat <- tmp.dat %>%
      mutate(day_diff = Date - lag(Date)) #for removing sites with missing days later
    message("Site ", nm, " has ", max(tmp.dat$day_diff, na.rm = TRUE), " day(s) at most between records.")
    
    write_csv(tmp.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))
    message("Flow dat for site ", nm, " downloaded.")
    
  }
}

rm(nm, tmp.dat)

## filter for 10 yr period with NO missing days----
missing_day_length <- 1 #set allowable amount of time between continuous data recordings (1 for no missing days)

#set function for saving info to dataframe
get_flow_dates <- function(dat) { 
  max_diff <- as.numeric(max(dat$day_diff, na.rm = TRUE))
  years <- as.numeric(dat$Date[nrow(dat)] - dat$Date[1]) / 365
  start_date <- dat$Date[1]
  end_date <- dat$Date[nrow(dat)]
  
  return(list(max_diff = max_diff, years = years, start_date = start_date, end_date = end_date))
}

#set function for clipping based on index location
clip_time <- function(dat, index_index, period_index) { 
  if(index_index == 1) {
    return(dat[1:(period_index[index_index]-1), ] %>% mutate(day_diff = Date - lag(Date)))
    
  } else if(index_index == length(period_index)) {
    return(dat[period_index[index_index-1]:period_index[index_index], ] %>% mutate(day_diff = Date - lag(Date)))
    
  } else {
    return(dat[period_index[index_index-1]:(period_index[index_index]-1), ] %>% mutate(day_diff = Date - lag(Date)))
    
  }
}

#set df for saving info on each site
site_check <- data.frame(site_no = site_nm, action = NA, 
                         l.max_diff = NA, l.years = NA, l.start_date = NA, l.end_date = NA, #longest period
                         n.max_diff = NA, n.years = NA, n.start_date = NA, n.end_date = NA) #newest period > 10 yrs

for(i in seq_len(nrow(site_check))) {
  tmp.dat <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check$site_no[i], ".csv"), col_types = cols(site_no = col_character()))
  
  if(max(tmp.dat$day_diff, na.rm = TRUE) <= missing_day_length) { #if whole dataset has no missing consecutive days, keep it all/no change
    site_check$action[i] <- "keep all"
    
    m <- get_flow_dates(tmp.dat)
    site_check[i, paste0("l.", names(m))] <- m
    
    message("Site ", site_check$site_no[i], " meets all time requirements, no change to data.")
    
  } else {
    
    ind <- which(tmp.dat$day_diff > missing_day_length) #get index of all time period jumps
    ind[length(ind)+1] <- nrow(tmp.dat) #add for getting last chunk
    
    t <- numeric(length(ind)) #prep
    
    #figure out # continuous years between > 1 day breaks
    for (j in seq_along(ind)) {
      if(j == 1) { #for first instance, find difference from start of records
        t[j] <- as.numeric(tmp.dat$Date[ind[j]-1] - tmp.dat$Date[1]) / 365
      } else {
        if(j == length(ind)) { #for the last one, include index row
          t[j] <- as.numeric(tmp.dat$Date[ind[j]] - tmp.dat$Date[ind[j-1]]) / 365
        } else {
          t[j] <- as.numeric(tmp.dat$Date[ind[j]-1] - tmp.dat$Date[ind[j-1]]) / 365
        }
      }
    }
    
    if(all(t < min_yr_period)) { #there is not enough data to pull from, remove from set
      site_check$action[i] <- "remove"
      message("Site ", site_check$site_no[i], " does not have ", min_yr_period, " yrs of continuous data.")
      next #go to next site
    }
    
    site_check$action[i] <- "clip"
    
    x <- which(t >= min_yr_period)
    
    long <- which(t == max(t[x])) #get longest period with > 10 yrs of data
    new <- tail(x, 1) #get most recent period with >10 yrs of data
    
    long.dat <- clip_time(tmp.dat, long, ind)
    m <- get_flow_dates(long.dat)
    site_check[i, paste0("l.", names(m))] <- m
    
    write_csv(long.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check$site_no[i], "_long.csv"))
    
    if(long != new) { #if there is a more recent period that isn't as long, save that too
      new.dat <- clip_time(tmp.dat, new, ind)
      m <- get_flow_dates(new.dat)
      site_check[i, paste0("n.", names(m))] <- m
      
      write_csv(new.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check$site_no[i], "_new.csv"))
    }
  
    message("Site ", site_check$site_no[i], " had > ", missing_day_length, " missing contig. days, long + new sections saved.")
    
  } #end else section for > 1 day breaks
}

site_check <- site_check %>% mutate(across(contains("date"), as.Date))
write_csv(site_check, paste0(PATH, "/02_EnvDat/raw_daily_flow/daily_flow_data_info.csv"))

rm(list = ls())

# HIT metric calc ----
PATH <- getwd()

#read in info on gages + flow data previously downloaded
site_check <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/daily_flow_data_info.csv"), col_types = cols(site_no = col_character())) %>%
  mutate(across(contains("date"), as.Date))

#which variables to pull with adjusted discharge (+ 0.01)
#redo for (ma4, ma9, ma10, ma11, mh18?, mh19?, fh11?, dh22?, dh23?, dh24?, tl3?, th3?, ra6, ra7) - last two automatic adj. in function
vars <- c("ma4", "ma9", "ma10", "ma11", "mh18", "mh19", "fh11", "dh22", "dh23", "dh24", "tl3", "th3")
# vars <- c("ma4", "ma9", "ma10", "ma11", "mh18", "mh19", "fh11", "dh22", "dh23", "dh24", "tl3", "tl4", "ta1", "ta2", "ta3", "th3")
#potentially just ma4, ra6, and ra7 overall though

hit_df <- data.frame()
# hit_df <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240122.csv"))

error_df <- data.frame(site_no = NA, error = NA)

for(i in seq_len(nrow(site_check))) {
  if(site_check$action[i] == "remove") { next } #skip those no longer being used bc of lack of data
  if(site_check$site_no[i] %in% hit_df$site_no) { next } #skip those already run for debugging reruns
  
  #read in flow dataset prev. downloaded
  site <- site_check$site_no[i]
  f <- file.path(PATH, "02_EnvDat/raw_daily_flow", 
                 paste0(site, ifelse(file.exists(file.path(PATH, "/02_EnvDat/raw_daily_flow", paste0(site, "_long.csv"))), #check if clipped file exists
                   "_long.csv", ".csv")))
  
  q <- read_csv(f, col_types = cols(site_no = col_character(), Date = col_date(), X_00060_00003 = col_number())) %>%
    mutate(month = format(Date, "%m"), day = format(Date, "%d"))
  
  if(any(grepl("OUTFLOW_00060", names(q)))) {     
    tmp <- data.frame(site_no = site, error = as.character("INFLOW/OUTFLOW column exists."), stringsAsFactors = FALSE) #make error for df
    error_df <<- bind_rows(error_df, tmp) #save error for later
    next 
    } #skip springs with inflow/outflow measurements
  
  #get complete water years (find first oct. 1, last sep. 30, clip everything after) + validate data
  ind1 <- head(which(q$month == "10" & q$day == "01"), 1)
  ind2 <- tail(which(q$month == "09" & q$day == "30"), 1)
  
  q.clean <- validate_data(as.data.frame(q[ind1:ind2, c("Date", "X_00060_00003")]), yearType = "water")
  
  #get drainage area + flood recurrence thresh
  site.info <- readNWISsite(siteNumber = site) #drainage area
  
  pk.fl <- readNWISpeak(siteNumbers = site, startDate = min(q.clean$date), endDate = max(q.clean$date))
  if(nrow(pk.fl) == 0) { #can't get peak
    pk.fl <- NULL
    fl.t <- NULL
  } else {
    pk.fl <- pk.fl %>% #peak flow
      filter(!is.na(peak_va) & !is.na(peak_dt)) #can't 'get_peakthreshold' if there are NA values in df
    
    if(any(pk.fl$peak_va == 0)) { #adjustment for log calc
       pk.fl$peak_va <- pk.fl$peak_va + 0.01 
       } 
    fl.t <- get_peakThreshold(q.clean[c("date", "discharge")], #flood recurrence threshold
                              pk.fl[c("peak_dt", "peak_va")])
  }
  
  ##calc HIT (all) ----
  tryCatch({
    hit_all <- calc_allHIT(q.clean, 
                           yearType = "water",
                           drainArea = site.info$drain_area_va[1],
                           floodThreshold = fl.t)
  }, error = function(e) {
    tmp <- data.frame(site_no = site, error = as.character(e), stringsAsFactors = FALSE) #get error
    error_df <<- bind_rows(error_df, tmp) #save error for later
  })
  
  ##calc HIT + 1 for log10 metrics ----
  q.clean2 <- q.clean %>% mutate(discharge = discharge + 0.01)
  if(!is.null(pk.fl)) {
    fl.t <- get_peakThreshold(q.clean2[c("date", "discharge")], #updated flood recurrence threshold
                              pk.fl[c("peak_dt", "peak_va")]) 
  }
  
  tryCatch({
    hit_all_2 <- calc_allHIT(q.clean2, 
                             yearType = "water",
                             drainArea = site.info$drain_area_va[1],
                             floodThreshold = fl.t)
  }, error = function(e) {
    tmp <- data.frame(site_no = site, error = paste0("adj: ", as.character(e)), stringsAsFactors = FALSE) #get error
    error_df <<- bind_rows(error_df, tmp) #save error for later
  })
  
  ix <- which(error_df$site_no == site)
  if(length(ix) >= 2) { next } #both output errors
  if(length(ix) == 1) {
    if(!grepl("adj", error_df$error[ix])) { #if only the adjusted run worked
      hit_all <- hit_all_2
    } #else is just hit_all by itself anyway
  }
  
  if(!site %in% error_df$site_no) {
    #add adjusted variables in
    hit_all_2 <- hit_all_2 %>%
      filter(indice %in% vars) %>% mutate(indice = paste0(indice, "_adj"))
    hit_all <- bind_rows(hit_all, hit_all_2)
  }

  #and mag7
  tmp <- calc_magnifSeven(q.clean, yearType = "water", digits = 3)
  
  #combine and add to full gage dataframe
  hit_all <- bind_rows(hit_all, tmp) %>%
    pivot_wider(names_from = indice, values_from = statistic) %>%
    mutate(site_no = site,
           hit_start_date = min(q.clean$date),
           hit_end_date = max(q.clean$date),
           hit_ttl_yrs = round(as.numeric(hit_end_date - hit_start_date) / 365, 3))
  
  hit_df <- bind_rows(hit_df, hit_all)
  
  write_csv(hit_df, paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_", gsub("-", "", Sys.Date()), ".csv"))
  message(round(i/nrow(site_check)*100, 2), "% complete")
}

error_df <- error_df %>% distinct() %>% filter(!is.na(site_no))
write_csv(error_df, paste0(PATH, "/02_EnvDat/raw_daily_flow/usgs_gage_hit_errors_", gsub("-", "", Sys.Date()), ".csv"))

hit_df <- hit_df %>%
  relocate(site_no, hit_start_date, hit_end_date, hit_ttl_yrs)
write_csv(hit_df, paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_", gsub("-", "", Sys.Date()), ".csv"))

rm(hit_all, pk.fl, q, q.clean, site.info, tmp, f, fl.t, ind1, ind2, site)

# get NHD stream code ----
#assign gages COMIDs
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv")) %>% #gage location
  filter(site_no %in% hit_df$site_no)
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp")) %>% #nhd streams prev. clipped to region (good to use clipped because of MEM usage)
  st_transform(5070)

#convert gage locations to spatial w/ projected crs
temp.sites <- st_as_sf(g.info %>% select(site_no, dec_lat_va, dec_long_va), coords = c("dec_long_va", "dec_lat_va"), crs = 4269, remove = F)
alb.sites <- temp.sites %>% st_transform(5070)

#assign COMID to each site + calc distance to nearest line (meters)
nr.line <- st_nearest_feature(alb.sites, nhd, check_crs = TRUE) #find nearest nhd strm
nhd_near <- nhd[nr.line,]
g.info$COMID <- nhd_near$COMID #add comid match to sites
dist <- as.vector(st_distance(alb.sites, nhd_near, by_element = TRUE)) #get distance to nearest nhd strm
g.info$dist2strm_m_nhd <- dist

#checking distance worked
# ggplot() +
#   geom_sf(data = nhd %>% filter(COMID == "22700332")) +
#   geom_sf(data = alb.sites %>% filter(site_no == "07360800")) +
#   ggspatial::annotation_scale() +
#   theme_minimal()

# get flow regime ----
#assign flow_type to each site + calc distance to nearest line (meters)
flw <- st_read(dsn = paste0(PATH, "/02_EnvDat/raw_flow_regime_dat/Interior Highlands Natural Flow Regimes.gdb"), layer = "Polylines") %>% select(-PopupInfo) %>%
  st_transform(5070)
flw <- flw %>% filter(Flow_type != "BigR")
nr.line <- st_nearest_feature(alb.sites, flw, check_crs = TRUE) #find nearest strm not river
flw_near <- flw[nr.line,]
g.info$flw_name <- flw_near$Name #add strm names
g.info$flw_type <- flw_near$Flow_type #add flow types
dist <- as.vector(st_distance(alb.sites, flw_near, by_element = TRUE)) #get distance from site to strm assigned
g.info$dist2strm_m_flw <- dist

rm(nhd, flw, alb.sites, flw_near, nhd_near, dist, nr.line, temp.sites)

#save location info
g.info <- g.info %>%
  left_join(., hit_df %>% select(site_no, contains("hit_"))) %>%
  mutate(across(contains("dist2strm"), ~ round(.x, digits = 3)))
write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))

# get HDI of gages ----
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
g.ii <- st_read(dsn = paste0(PATH, "/02_EnvDat/GAGES_II/point_shapefile"), layer = "gagesII_9322_sept30_2011") %>%
  select(STAID, DRAIN_SQKM) %>%
  rename(drain_sqkm = DRAIN_SQKM) %>%
  st_drop_geometry()
g.ii.hdi <- read_csv(paste0(PATH, "/02_EnvDat/GAGES_II/basinchar_and_report_sept_2011/spreadsheets_in_csv_format/conterm_bas_classif.txt")) %>%
  select(STAID, CLASS, HYDRO_DISTURB_INDX) %>%
  rename(gage_class = CLASS, HDI = HYDRO_DISTURB_INDX)

g.info <- g.info %>%
  left_join(., g.ii, by = c("site_no" = "STAID")) %>%
  left_join(., g.ii.hdi, by = c("site_no" = "STAID")) %>%
  select(-c(agency_cd, srs_id, access_cd, medium_grp_cd, parm_grp_cd, loc_web_ds)) 
write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))

g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))


#check + fix errors ----
#ignore quantile error and inflow/outflow skips
#freq high throws the error, duration high also, potentially has to do with the all flow 0 values in 1990
#timing low, timing high
#for all, if flood threshold is removed, values are returned, same for duration high, though some will be NA
# calc_rateChange(q.clean, 
#             yearType = "water",
#             drainArea = site.info$drain_area_va[1],
#             floodThreshold = fl.t
#             )

#compare to old HIT ----
# old.hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv"),
#                     col_types = cols(STAID_0 = col_character())) %>%
#   select(-c(FID, STAID...2, STANAME, `Years in Op`, `Reference?`, `Start Date`, `End Date`, STAID...13), -contains("GAGESII"), -contains("Active"))
# 
# tmp <- left_join(hit_df, old.hit, by = c("site_no"="STAID_0")) %>%
#   # left_join(old.hit, hit_df, by = c("STAID_0"="site_no")) %>% rename(site_no = STAID_0) %>%
#   relocate(site_no, contains("date"), contains("year"), matches(paste(names(hit_df)[-c(1:3)], collapse = "|"))) %>%
#   mutate(across(ends_with(".x"), ~ . - get(sub("\\.x$", ".y", cur_column())), .names = "{col}_diff")) %>%
#   select(site_no, contains("date"), contains("year"), sort(names(.)))

#     range standardize hydrologic data for gages (by dividing each metric by their absolute values to transform metrics to comparable scales?)

