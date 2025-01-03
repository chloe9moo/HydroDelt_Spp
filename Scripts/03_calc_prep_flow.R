## FLOW CALC + ESTIMATION ##
## doing this before filtering occurrences b/c need date range of flow data

library(tidyverse); library(sf)
library(dataRetrieval); library(EflowStats)
options(readr.show_col_types = FALSE)

PATH <- getwd()

# pull daily flow ----
min_yr_period <- 15 #set number of years req. for use in analyses

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

#latest date for flow data set
occ <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_inthigh_long_20240716.csv"))
max.date <- max(occ$date, na.rm = TRUE)
rm(occ)

#pull flow data and save
for(i in seq_along(site_nm)) { 
  nm <- site_nm[[i]]
  
  if(file.exists(paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))) { #in case it fails at any point
    message("Flow dat for ", nm, " already exists."); next
  } else {
    
    #set last pull date
    end.date <- max(flow.gage[flow.gage$site_no == nm, ]$end_date)
    
    if(end.date > max.date) { #if later than occurrence data, then pull up to that date
      end.date <- max.date
    }
    
    tmp.dat <- readNWISdv(siteNumbers = nm,
                          parameterCd = "00060",
                          statCd = "00003",
                          startDate = min(flow.gage[flow.gage$site_no == nm, ]$begin_date),
                          endDate = end.date)
    
    tmp.dat <- tmp.dat %>%
      mutate(day_diff = Date - lag(Date)) #for removing sites with missing days later
    message("Site ", nm, " has ", max(tmp.dat$day_diff, na.rm = TRUE), " day(s) at most between records.")
    
    write_csv(tmp.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))
    message("Flow dat for site ", nm, " downloaded.")
    
  }
}

rm(nm, tmp.dat)

## filter for 15 yr period with NO missing days----
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
                         n.max_diff = NA, n.years = NA, n.start_date = NA, n.end_date = NA, #newest period > 10 yrs
                         neg_val_ct = NA) #number of negative discharge values in the data 

for(i in seq_len(nrow(site_check))) {
  tmp.dat <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check$site_no[i], ".csv"), col_types = cols(site_no = col_character()))
  
  #handle negative discharge values (turn into NAs and then clip around them)
  if(any(names(tmp.dat) %in% c("X_00060_00003"))) { #don't worry about those INFLOW/OUTFLOW cols, bc not going to use those gages
    
    site_check[i, ]$neg_val_ct <- length(which(tmp.dat$X_00060_00003 < 0))
    if(any(tmp.dat$X_00060_00003 < 0)) {
      tmp.dat <- tmp.dat %>%
        mutate(X_00060_00003 = case_when(X_00060_00003 < 0 ~ NA,
                                         T ~ X_00060_00003)) %>%
        filter(!is.na(X_00060_00003)) %>%
        mutate(day_diff = Date - lag(Date))
    }
    
  }
  
  #if whole dataset has no missing consecutive days, keep it all/no change
  if(max(tmp.dat$day_diff, na.rm = TRUE) <= missing_day_length) { 
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
    
    long <- which(t == max(t[x])) #get longest period with > X yrs of data
    new <- tail(x, 1) #get most recent period with > X yrs of data
    
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

site_check <- site_check %>% 
  mutate(across(contains("date"), as.Date),
         action = ifelse(!is.na(l.years) & l.years < min_yr_period, "remove", action))

#used to check changes made on 1JUL2024, no major differences besides for those few sites with neg to NA values
# tmp_check <- left_join(daily_flow_data_info, site_check, by = "site_no")
# tmp_check <- tmp_check %>% 
#   mutate(across(ends_with(".x") & is.numeric & -contains("change"), ~ . - get(sub("\\.x", "\\.y", cur_column())), .names = "change_{col}")) %>% 
#   relocate(matches(names(site_check)))

write_csv(site_check, paste0(PATH, "/02_EnvDat/raw_daily_flow/daily_flow_data_info.csv"))

rm(list = ls())

# HIT metric calc ----
PATH <- getwd()

#read in info on gages + flow data previously downloaded
site_check <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/daily_flow_data_info.csv"), col_types = cols(site_no = col_character())) %>%
  mutate(across(contains("date"), as.Date))

#which variables to also pull with adjusted discharge (+ 0.01)
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
  #check and set which clipped file to load (if exists)
  if(file.exists(file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, "_long.csv")))) {
    f <- file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, "_long.csv"))
  } else {
    f <- file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, ".csv"))
  }
  #if newer clipped file exists, overwrite
  if(file.exists(file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, "_new.csv")))) { 
    f <- file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, "_new.csv"))
  }
  
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
  
  ##calc HIT + 0.01 for log10 metrics ----
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
  
  ## adjust some variable calculations ----
  #adjust dl17 (variability in low pulse duration) - returns NA when no events below low pulse occur (high 0 flow sites essentially)
  #replace with 0 to indicate no low pulse events
  if(any(is.na(hit_all[hit_all$indice == "dl17", ]$statistic))) {
    hit_all[hit_all$indice == "dl17", ]$statistic <- 0
  }
  
  #adjust ma24-35 (mean monthly CV) to ignore NA vals, rather than returning an NA value
  #the MHIT (MATLAB version of the package) did it this way, modified the source code for calc_magAverage
  if(any(q.clean$discharge == 0)) {
    ma24.35 <- q.clean %>%
      mutate(month_val = month(date)) %>%
      group_by(year_val, month_val) %>%
      summarise(meanFlow = mean.default(discharge),
                cvMonth = sd(discharge)/meanFlow) %>%
      group_by(month_val) %>%
      summarise(meanCV = mean.default(cvMonth, na.rm = TRUE)) #<< mod here
    ma24.35 <- round(ma24.35$meanCV * 100, digits = 3)
    
    ma24.35 <- data.frame(indice = c(paste0("ma", 24:35)),
                          statistic = ma24.35,
                          stringsAsFactors = FALSE)
    
    hit_all <- bind_rows(hit_all[!hit_all$indice %in% ma24.35$indice, ], ma24.35)
    
    rm(ma24.35)
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
  rm(hit_all, hit_all_2)
  message(round(i/nrow(site_check)*100, 2), "% complete")
}

error_df <- error_df %>% distinct() %>% filter(!is.na(site_no))
write_csv(error_df, paste0(PATH, "/02_EnvDat/raw_daily_flow/usgs_gage_hit_errors_", gsub("-", "", Sys.Date()), ".csv"))

#missing sites
# missing_site <- site_check[!site_check$site_no %in% hit_df$site_no, ]
# inf.check <- as.data.frame(which(is.infinite(as.matrix(hit_df[, sapply(hit_df, is.numeric)])), arr.ind = TRUE))

#prep + save
hit_df <- hit_df %>%
  relocate(site_no, hit_start_date, hit_end_date, hit_ttl_yrs) %>%
  mutate(across(everything(), ~ case_when(is.infinite(.x) ~ NA,
                                          is.nan(.x) ~ NA,
                                          T ~ .x)))
write_csv(hit_df, paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_", gsub("-", "", Sys.Date()), ".csv"))

rm(pk.fl, q, q.clean, site.info, tmp, f, fl.t, ind1, ind2, site, q.clean2)

###summarize missing data + outlier abundance ----
#pull out hit vars only
df <- hit_df[, sapply(hit_df, is.numeric) & !names(hit_df) %in% "hit_ttl_yrs"]
#NAs by site
na_count <- data.frame(name = hit_df$site_no, type = "site", na_count = NA, na_pct = NA)
na_count$na_count <- apply(df, 1, function(x) sum(is.na(x)))
na_count$na_pct <- round((na_count$na_count / ncol(df)) * 100, digits = 2)

#flag outliers by site
is_outlier <- function(x) {
  stdev <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  return(abs(x - m) > 3 * stdev) #greater than 3 standard deviations = outlier
}
outlier_id <- df %>%
  mutate(across(everything(), ~ is_outlier(.x), .names = "outflag_{col}"),
         outlier_total = rowSums(across(starts_with("outflag")), na.rm = TRUE),
         outlier_pct = round(outlier_total / ncol(df) * 100, digits = 2))

write_csv(outlier_id, paste0(PATH, "/02_EnvDat/raw_daily_flow/HIT_outlier_summary_20240813.csv"))

#add tracking columns to na_count
na_count$outlier_total <- outlier_id$outlier_total
na_count$outlier_pct <- outlier_id$outlier_pct

#plotting to check things
# ggplot() +
#   # geom_density(data = df, aes(ma20))
#   # geom_point(data = df %>% filter(!is.na(ma4)), aes(x = ma4, y = ma4_adj), size = 4)
#   geom_boxplot(data = outlier_id, aes(x = "ma4", y = ma4)) +
#   geom_point(data = outlier_id, aes(x = "ma4", y = ma4, color = outflag_ma4))

#NAs by variable
tmp <- data.frame(name = names(df), type = "HIT", na_count = NA, na_pct = NA)
tmp$na_count <- apply(df, 2, function(x) sum(is.na(x)))
tmp$na_pct <- round((tmp$na_count / nrow(df)) * 100, digits = 2)

#outliers across sites within variable
o <- outlier_id %>%
  select(contains("outflag")) %>%
  pivot_longer(everything(), names_to = "name", values_to = "outlier") %>%
  group_by(name) %>%
  summarise(outlier_total = sum(outlier, na.rm = TRUE)) %>%
  mutate(name = gsub("outflag_", "", name),
         outlier_pct = round(outlier_total / nrow(df) * 100, digits = 2))
tmp <- left_join(tmp, o)

na_count <- bind_rows(na_count, tmp)

write_csv(na_count, paste0(PATH, "/02_EnvDat/raw_daily_flow/HIT_missing_outlier_summary.csv"))

rm(df, o, outlier_id, tmp)

# get NHD stream code + flow regime ----
##NOTE: if redoing the HIT calc, don't need to do this section and HDI section bc it is attached to gage by site_no, not related to HIT calc
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv")) %>% #gage location
  filter(site_no %in% hit_df$site_no)

source(paste0(PATH, "/Scripts/XX_strm_flw_info_func.R"))

g.info <- get_strm_flow_info(g.info[, names(g.info) %in% c("site_no", "dec_long_va", "dec_lat_va")], 
                             id_col = c("site_no"), #for saving site ID
                             coord_cols = c("dec_long_va", "dec_lat_va"), #in correct order  (x, y)
                             site_crs = 4269, #crs of the original coordinate data
                             get_COMID = TRUE,
                             get_COMID_dist = TRUE,
                             get_flw_type = TRUE,
                             get_flw_dist = TRUE)

#save location info
g.info <- g.info %>%
  left_join(., hit_df %>% select(site_no, contains("hit_"))) %>%
  left_join(., na_count %>% filter(type == "site") %>% select(-type), by = c("site_no" = "name")) %>%
  mutate(across(contains("dist2strm"), ~ round(.x, digits = 3))) %>%
  rename(hit_na_count = na_count, hit_na_pct = na_pct)
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
  left_join(., g.ii.hdi, by = c("site_no" = "STAID"))
write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))

rm(list = ls())

#decide which gages // HIT vars to use ----
PATH <- getwd()

g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))

#which gages to keep?
g.info <- g.info %>%
  mutate(final_filter = case_when(hit_ttl_yrs < 15 ~ "remove", #some sites after clipping to meet continuous time req. have less than 10 yrs of data now
                                  dist2strm_m_nhd > 1000 ~ "remove", #if gage is far from NHD stream (> 1 km), (which is used for ID to many env vars), remove
                                  dist2strm_m_flw > 1000 ~ "remove", #if gage is far from classified flow regimes (> 1 km), remove
                                  #sites with missing data in target HIT variables
                                  T ~ "keep"))

table(g.info$flw_type) #evenish distribution of flow types?

#map to look at where gages are bein removed
#background shapes first
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop")
armook <- st_as_sf(maps::map("state", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>%
  st_transform(., crs = 4269)

sf.g <- g.info %>% st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269) %>% #EPSG 4269 = NAD83
  mutate(final_filter = factor(final_filter, levels = c("keep", "remove"))) %>%
  arrange(final_filter)
ggplot() +
  geom_sf(data = armook, fill = "lightgray", color = "black") +
  geom_sf(data = eco, fill = "darkgray", color = "black") +
  geom_sf(data = sf.g, aes(fill = final_filter, shape = flw_type, size = final_filter)) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_size_manual(values = c(3, 5)) +
  theme_bw()

#save for filtering later
write_csv(g.info, paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))

#data for making decisions 
out_info <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/HIT_outlier_summary.csv"))
hit <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240702.csv"))
var_table <- read_csv(paste0(PATH, "/02_EnvDat/environmental_variable_info.csv"))
na_count <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/HIT_missing_outlier_summary.csv"))
var_table <- left_join(var_table, na_count, by = c("variable_col_name" = "name"))
var_table <- var_table[!is.na(var_table$type),]

#which HIT variables to keep?
out_info$remove_gage <- g.info$final_filter
df <- out_info %>%
  filter(remove_gage == "keep") %>%
  select(matches(paste0("^", var_table[!is.na(var_table$include_hydro_fm24), ]$variable_col_name, "$")),
         matches(paste0("^", var_table[!is.na(var_table$include_hydro_fm24), ]$variable_col_name, "_adj$")))
tmp <- data.frame(name = names(df), type = "HIT", na_count = NA, na_pct = NA)
tmp$na_count <- apply(df, 2, function(x) sum(is.na(x)))
tmp$na_pct <- round((tmp$na_count / nrow(df)) * 100, digits = 2)

out_info %>%
  mutate(fm24_na_ct = pmap_int(across(
    matches(paste0("^", var_table[!is.na(var_table$include_hydro_fm24), ]$variable_col_name, "$"))
  ), ~ sum(is.na(c(...)))),
  fm24_na_pct = fm24_na_ct / length(var_table[!is.na(var_table$include_hydro_fm24), ]$variable_col_name) * 100,
  fm24_na_pct = round(fm24_na_pct, digits = 2)) %>%
  View()

View(out_info %>% select(ma4, ma4_adj))

#NOTE from 15JUN2024: removing ma4 from variable list, including adjusted ma4 - too many questionable values on edges (e.g., -400000 ???), made edit to 
##env variable info table under column include_hydro_fm24

#notes on errors from HIT calc----
#ignore quantile error and inflow/outflow skips
#ignore error in 1:sum(runLengths$values == T): NA/NaN argument
#freq high throws the error, duration high also, potentially has to do with the all flow 0 values in 1990
#timing low, timing high
#for all, if flood threshold is removed, values are returned, same for duration high, though some will be NA
# calc_rateChange(q.clean, 
#             yearType = "water",
#             drainArea = site.info$drain_area_va[1],
#             floodThreshold = fl.t
#             )

#compare to old HIT ----
# g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv"))
# hit_df <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240715.csv"))
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

