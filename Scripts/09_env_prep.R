## OBTAIN + PREP ENV VARIABLES FOR GF MODELS ##

library(tidyverse); library(sf); library(terra); library(exactextractr); library(parallel)
options(readr.show_col_types = FALSE, dplyr.summarise.inform=F)

PATH <- getwd()

#load in sites for attaching env var data
occ.sites <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
occ.sites <- lapply(occ.sites, read_csv, col_types = cols(lat = col_number(),
                                                          long = col_number(),
                                                          .default = "c"))

#get COMIDs (needed to obtain some env vars)
comids <- lapply(occ.sites, function(x) select(x, COMID))
comids <- bind_rows(comids)
comids <- unique(comids$COMID)

# StreamCat Vars ----
library(StreamCatTools)

#select variables to pull
# sc.names <- sc_get_params(param = "name") #variable names
# 
# sc.name.df <- data.frame()
# for(metric in sc.names) { #find full names for matching
#   fullname <- sc_fullname(metric = metric)
#   
#   if(length(fullname) == 0) { fullname <- NA }
#   
#   fn1 <- data.frame(env_name = metric, full_name = fullname)
#   
#   sc.name.df <- bind_rows(sc.name.df, fn1)
# }
# rm(fn1, fullname, metric)
# write_csv(sc.name.df, paste0(PATH, "/02_EnvDat/StreamCat/stream_cat_variable_names.csv"))

##added column manually to select vars, based on previous GF work
sc.name.df <- read_csv(paste0(PATH, "/02_EnvDat/StreamCat/stream_cat_variable_names.csv"))
sc.name.df <- sc.name.df[sc.name.df$include == "yes" & !is.na(sc.name.df$include), ]

#pull data
metric.v <- sc.name.df$env_name
metric.v <- paste0(metric.v, collapse = ",") #format for function

comids.l <- split(comids, ceiling(seq_along(comids) / 50)) #split vector of comids bc it's too long for a single pull

stream_cat <- vector("list", length = length(comids.l))
for(i in seq_along(comids.l)) { 
  
  stream_cat[[i]] <- sc_get_data(metric = metric.v, aoi = "watershed,other", comid = comids.l[[i]])
  
  message(round(i/length(comids.l)*100, digits = 2), "% complete")
  Sys.sleep(1)
  
  }

stream_cat <- bind_rows(stream_cat)

m.c <- comids[!comids %in% stream_cat$COMID] #comids not returned

stream_cat_sum <- stream_cat %>%
  pivot_longer(-COMID, names_to = "var", values_to = "value") %>%
  group_by(var) %>%
  summarise(mn = mean(value, na.rm = T),
            med = median(value, na.rm = T),
            max = max(value, na.rm = T),
            min = min(value, na.rm = T),
            cv = (sd(value, na.rm = T) / mean(value, na.rm = T)),
            prop_na = sum(is.na(value))/nrow(stream_cat)*100) %>%
  mutate(across(-var, ~ round(.x, digits = 3)))

#save
write_csv(stream_cat, file = paste0(PATH, "/02_EnvDat/StreamCat/stream_cat_vars_all_sites.csv"))
write_csv(stream_cat_sum, file = paste0(PATH, "/02_EnvDat/StreamCat/stream_cat_var_summary.csv"))

rm(sc.name.df, stream_cat, stream_cat_sum, i, m.c, metric.v)

# Env Change Over Time variables ----
#obs - exp / exp, where exp is 'normal' or average across full set of years (??), obs is monthly/seasonal avg, then average those values (??)
#sum of months(obs - exp / exp) = seasonality index per McManamay et al. 2022 (hydro alt dataset)

## prep ----
#get date range of occurrence data
# load in dat with dates (if available)
# file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "inthigh_long", full.names = TRUE)
# occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
#                                                          long = col_number(),
#                                                          lat_new = col_number(),
#                                                          long_new = col_number(),
#                                                          taxa_count = col_number(),
#                                                          .default = "c"))
# 
# lapply(occ.list, function(x) {
#   #filter to match wide data parameters
#   x <- x %>%
#     filter(dist2gage_m_15yr <= 15000 & dist2strm_m_flw <= 5000) %>% #make sure within 15km of a gage + 5 km classified flow line
#     filter(!gage_no_15yr %in% c("06906800", "07020550", "07332500"))
#   #combine all dates
#   all.dates <- c(x$date_min, x$date_max, x$date)
#   all.dates <- all.dates[!is.na(all.dates)]
# 
#   paste0(min(all.dates), " to ", max(all.dates))
#   # return(all.dates)
# })

# fish.dates <- all.dates[[2]]
# fish.dates <- as.Date(fish.dates)
#based on looking into the data a bit more, the ADEQ sites with a date of 1900 are likely wrong (looked at ADEQ website)

#bugs: Aug. 1964 to Jan. 2022
#fish: Jan. 1923 to Jun. 2021

#get watersheds
huc <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/huc12_interior_highlands_crop.shp"))

#pair watersheds w/ comids
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp"))
nhd <- nhd[nhd$COMID %in% comids,]

if(file.exists(paste0(PATH, "/02_EnvDat/HUCs_NHDs/nhd_huc_intersection_info.csv"))) { #load intersection info if starting over
  
  huc_nhd_df <- read_csv(paste0(PATH, "/02_EnvDat/HUCs_NHDs/nhd_huc_intersection_info.csv"))
  
} else { #get intersection info if starting over
  
  huc_nhd_int <- st_intersects(huc, nhd) #returns list where each numbered element = huc feature in order, in element is feature number intersecting from nhd
  
  huc_nhd_df <- data.frame(comid_idx = seq(1, nrow(nhd)), comid = NA, huc_id = NA, prop_in_huc = NA)
  
  for(i in seq_len(nrow(nhd))) {
    huc_nhd_df[i,]$comid <- nhd[i,]$COMID #get nhd index COMID
    
    int_idx <- sapply(huc_nhd_int, function(x) { i %in% x }) #find where the intersection is
    
    if(length(which(int_idx)) > 1) { #if the stream overlaps more than one watershed..
      nhd_length <- st_length(nhd[i, ])
      
      prop_length <- vector("numeric", length(which(int_idx)))
      
      for(j in seq_along(which(int_idx))) {
        nhd_inter <- st_intersection(nhd[i,], huc[which(int_idx)[j], ]) #get side of linestring in huc
        
        nhd_inter <- st_length(nhd_inter) #get length in that huc
        
        prop_length[j] <- round(nhd_inter/nhd_length*100, digits = 3) #get proportion of stream in huc
      }
      
      huc_nhd_df[i, ]$prop_in_huc <- prop_length[which.max(prop_length)] #save proportion overlapped for later
      
      if(prop_length[which.max(prop_length)] < 70) { #i think it makes sense to pull date for more than 1 huc for some of these streams
        
        huc_nhd_df[i, ]$huc_id <- paste0(huc[which(int_idx), ]$huc12, collapse = "|")
        
      } else {
        
        huc_nhd_df[i, ]$huc_id <- huc[which(int_idx)[which.max(prop_length)], ]$huc12 #now get huc id for most overlapping watershed
        
      }
      
    } else {
      huc_nhd_df[i, ]$prop_in_huc <- 100 #if only one ws returned, assumed 100% strm in it
      
      huc_nhd_df[i, ]$huc_id <- huc[which(int_idx), ]$huc12
    }
    
    message(round(i/nrow(huc_nhd_df)*100, digits = 2), "% complete")
  }
  
  write_csv(huc_nhd_df, paste0(PATH, "/02_EnvDat/HUCs_NHDs/nhd_huc_intersection_info.csv"))
  
  #plot to check
  # ggplot() + geom_sf(data = huc[c(58, 48),]) + geom_sf(data = huc[58,], color = "red") + geom_sf(data = nhd[c(7),])
  
  rm(int_idx, nhd_length, prop_length, nhd_inter, huc_nhd_int, i, j)
  
}

# use nhd/huc overlap to find average for each month x year needed within watershed boundaries
huc_ids <- paste0(huc_nhd_df$huc_id, collapse = "|")
huc_ids <- strsplit(huc_ids, "\\|")[[1]]
huc_ids <- unique(huc_ids)

huc <- huc[huc$huc12 %in% huc_ids, ]

## air temp change ----
##note: prism data across almost all years were previously downloaded and cropped in 07_get_fstrm_temp_data.R
#get average temp within each watershed for each year of interest (Jan 1923 to Jan 2022)
prism_files <- list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters"), full.names = TRUE)

a_temp <- data.frame()
a_temp <- mclapply(prism_files, mc.cores = 8, function(x) { #~ 18 GB MEM at max using 8 cores
  date <-  gsub(".*?(\\d{6}).*", "\\1", x)
  
  if(nchar(date) != 6) { return(NULL) }
  
  pr1 <- rast(x) #load in clipped raster
  
  #extract temperature values at all gage sites
  a_temp1 <- exact_extract(pr1, huc, fun = 'weighted_mean', weights = 'area', force_df = TRUE, progress = FALSE) #exactextractr::exact_extract is faster than terra::extract
  
  names(a_temp1) <- "air_temp"
  #add hucID and date to df
  a_temp1$year <- gsub("^(.{4}).*", "\\1", date)
  a_temp1$month <- gsub(".*(.{2})$", "\\1", date)
  a_temp1 <- bind_cols(st_drop_geometry(huc[, "huc12"]), a_temp1)
  
  return(a_temp1)
})

a_temp <- bind_rows(a_temp)
a_temp <- a_temp %>% arrange(huc12, year, month)

write_csv(a_temp, paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_huc_lvl.csv"))

#calculate change
a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_huc_lvl.csv"))

#bugs: Aug. 1964 to Jan. 2022
#fish: Jan. 1923 to Jun. 2021

##annual deficit
calc_var_diff <- function(x, #variable dataframe that has, at minimum, env variable + time
                          var2grp = c("huc12", "year", "season"), #group levels, right now this only works across years, not within a single year
                          norm_comp = TRUE, #calc deficit from normal, departure from normal, percent of normal (all averaged over time)
                          pct_change = TRUE, #average percent change across time for each lowest group
                          delta_var = TRUE, #difference in var from start to end of period, grouped by N years
                          num_yr_start_end = 5, #for difference from start to end of period, how many years to pull from
                          c_var_diff = TRUE, #cumulative change over length of time (essentially length of the trend line)
                          growth_rate = TRUE #growth rate calculated from decomposed time series, not grouped by season or anything
                          ) {
  # x <- air_temp #to test
  all_res <- x[, names(x) %in% var2grp[!grepl("year", var2grp)]]
  all_res <- distinct(all_res)
  
  if(norm_comp) {
    #get normal = average value over all years ----
    norm <- x %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm = mean(air_temp))
    
    obs <- x %>%
      group_by(across(all_of(var2grp))) %>% 
      summarise(obs = mean(air_temp)) %>%
      ungroup() %>%
      left_join(., norm, by = names(.)[names(.) %in% names(norm)])
    
    #within years deficits ----
    def <- obs %>%
      mutate(def = (obs - norm) / norm) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_def = mean(def)) %>%
      ungroup()
    
    #departure from normal (obs - norm) ----
    depart <- obs %>%
      mutate(depart = obs - norm) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_depart = mean(depart)) %>%
      ungroup()
    
    #percentage of normal (obs/normal * 100) ----
    pct <- obs %>%
      mutate(pct = obs/norm*100) %>% 
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_pct = mean(pct)) %>%
      ungroup()
    
    all_res <- all_res %>%
      left_join(., def, by = names(.)[names(.) %in% names(def)]) %>%
      left_join(., depart, by = names(.)[names(.) %in% names(depart)]) %>%
      left_join(., pct, by = names(.)[names(.) %in% names(pct)])
  }
  
  
  #avg percent change ----
  if(pct_change) {
    mn.pch <- x %>% 
      group_by(across(all_of(var2grp))) %>%
      summarise(mn_temp = mean(air_temp)) %>%
      arrange(across(all_of(var2grp))) %>% 
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(pct_ch = (mn_temp / lag(mn_temp) - 1) * 100,
             pct_ch = ifelse(mn_temp > lag(mn_temp), abs(pct_ch), pct_ch)) %>% #second line to handle transitions from neg. to pos.
      summarise(mn_pch = mean(pct_ch, na.rm = TRUE))
    
    all_res <- all_res %>% left_join(., mn.pch, by = names(.)[names(.) %in% names(mn.pch)])
  }

  
  #diff between start and end of sample period ----
  if(delta_var) {
    dates <- unique(x$year)
    dates <- sort(dates)
    
    d_temp <- x %>%
      mutate(time_point = case_when(year %in% head(dates, num_yr_start_end) ~ "first",
                                    year %in% tail(dates, num_yr_start_end) ~ "last")) %>%
      filter(!is.na(time_point)) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)])), time_point) %>%
      summarise(avg_temp = mean(air_temp, na.rm = TRUE)) %>%
      arrange(across(all_of(var2grp[!grepl("year", var2grp)])), time_point) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(delta_temp = avg_temp - lag(avg_temp)) %>%
      select(any_of(var2grp), delta_temp) %>%
      filter(!is.na(delta_temp))
    
    all_res <- all_res %>% left_join(., d_temp,  by = names(.)[names(.) %in% names(d_temp)])
  }
  
  
  #cumulative diff. of over time ----
  if(c_var_diff) {
    c_temp <- x %>%
      group_by(across(all_of(var2grp))) %>%
      summarise(ssn_temp = mean(air_temp)) %>% #one seasonal value per year
      arrange(across(all_of(var2grp))) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(cdelt_temp = ssn_temp - lag(ssn_temp)) %>%
      summarise(cdelt_temp = sum(abs(cdelt_temp), na.rm = TRUE))
    
    all_res <- all_res %>% left_join(., c_temp,  by = names(.)[names(.) %in% names(c_temp)])
  }
  
  #ratio of delta var : cum. var ----
  if(all(c(delta_var, c_var_diff))) {
    all_res <- all_res %>%
      mutate(change_ratio = cdelt_temp / delta_temp)
  }
  
  #growth rate from decomposed time series ----
  if(growth_rate) {
    ts.l <- split(x, x$huc12)
    
    ts.l <- lapply(ts.l, function(y) {
      y <- y %>% arrange(date)
      
      ts.gr <- ts(y$air_temp, start = c(y[1, ]$year, y[1, ]$month), end = c(y[nrow(y), ]$year, y[nrow(y), ]$month), frequency = 12)
      
      m <- decompose(ts.gr, type = "multiplicative")
      trend <- na.omit(m$trend)
      trend <- diff(log(trend))
      
      #avg percent growth rate
      pct.avg.gr <- mean(trend) * 100
      
      res <- data.frame(huc12 = unique(y$huc12), pct_avg_gr = pct.avg.gr)
      
      return(res)
    })
    
    ts.l <- bind_rows(ts.l)
    
    all_res <- all_res %>% left_join(., ts.l,  by = names(.)[names(.) %in% names(ts.l)])
  }

  
  return(all_res)
}

summarize_temp_delt <- function(air_temp = a_temp, #temperature df just created, with HUC IDs and year x month
                                date_min = as.Date("1964-08-01"), #set taxa date min and max for filtering air temperature
                                date_max = as.Date("2022-01-31"),
                                seasonal = TRUE
) {
  
  #filter temperature to desired date range
  air_temp <- air_temp %>%
    mutate(date = as.Date(paste(year, month, "01", sep = "-")),
           season = case_when(month %in% c("12", "01", "02") ~ "winter",
                              month %in% c("03", "04", "05") ~ "spring",
                              month %in% c("06", "07", "08") ~ "summer",
                              month %in% c("09", "10", "11") ~ "fall")) %>%
    filter(date >= date_min & date <= date_max)
  
  #seasonal differences
  if(seasonal){
    output <- calc_var_diff(air_temp, var2grp = c("huc12", "year", "season"),
                            norm_comp = TRUE, pct_change = TRUE, 
                            delta_var = TRUE, num_yr_start_end = 5,
                            c_var_diff = TRUE, growth_rate = TRUE)
  }

  ssn <- pivot_wider(ssn, names_from = season, values_from = c(norm_def, norm_depart, norm_pct))
  
  TO DO:
    CHECK FUNCTION WORKS, FROM FRESH ENVIRONMENT (seems to work fine)
    RUN FOR BOTH TAXA, SAVE THE RESULTS
    UPDATE CODE TO WORK WITH OTHER VARIABLES (E.G. DO PRECIP NEXT)

}

#comparing diff change ratios
air_temp %>%
  #order is high neg, high pos, middle value
  filter((huc12 == "071401020407" & season == "winter") | (huc12 == "102901070301" & season == "winter") | (huc12 == "111402010501" & season == "winter")) %>%
  group_by(huc12) %>%
  arrange(date) %>%
  mutate(scaled_var = (air_temp - mean(air_temp)) / sd(air_temp)) %>%
  ggplot(aes(date, scaled_var)) +
  geom_line(aes(color = huc12), alpha = 0.3) +
  geom_smooth(method = "loess", span = 0.7, aes(color = huc12)) + #change span to edit smoothing amount (higher is more)
  theme_classic()



#from NOAA - https://water.weather.gov/precip/ (some ideas here)
#precip deficit (???)
#departure from normal? (not sure how this would work if we need 1 value for each site)
#average departure from normal? could also do amount of variation??


#stream temp change

#precip change



#land cover change

##ADD HDI VARIABLES FROM FOX & MAGOULICK 2022?







t <- nhd %>%
  mutate(missing_c = ifelse(COMID %in% m.c, "wrong", "right")) %>%
  group_by(missing_c) %>%
  summarise(across(where(is.numeric), mean))




