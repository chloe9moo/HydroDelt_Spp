## OBTAIN + PREP ENV VARIABLES FOR GF MODELS ##

library(tidyverse); library(sf); library(terra); library(exactextractr); library(parallel)
options(readr.show_col_types = FALSE, dplyr.summarise.inform=F)

PATH <- getwd()

#set which sections to run ----
get_stream_cat <- FALSE; get_hydro_alt <- FALSE; get_adj_hit <- FALSE
get_occ_date_range <- FALSE
get_huc_lvl_air_temp <- FALSE; summarize_air_temp <- TRUE
download_prism_ppt <- FALSE; clip_prism_ppt <- FALSE; get_huc_lvl_ppt <- FALSE; summarize_ppt <- TRUE
summarize_stream_temp <- TRUE
download_nlcd <- FALSE; get_huc_lvl_nlcd <- FALSE; summarize_nlcd <- TRUE

#load in sites for summarizing env var data
occ.sites <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
occ.sites <- lapply(occ.sites, read_csv, col_types = cols(lat = col_number(),
                                                          long = col_number(),
                                                          .default = "c"))

#get gage nos
gage_ids <- lapply(occ.sites, function(x) select(x, gage_no_15yr))
gage_ids <- bind_rows(gage_ids)
gage_ids <- unique(gage_ids$gage_no_15yr)

#get COMIDs (needed to obtain some env vars)
comids <- lapply(occ.sites, function(x) select(x, COMID))
comids <- bind_rows(comids)
comids <- unique(comids$COMID)

# StreamCat Vars ----
if(get_stream_cat) {
  
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

}

# Hydro Alt dataset ----
if(get_hydro_alt) {
  
#(McManamay et al. 2022)
h.alt <- read_csv(paste0(PATH, "/02_EnvDat/hydro_alt_disturb/predicted_alteration_US_model.csv"), col_types = cols(COMID = col_character(),
                                                                                                                   .default = "n"))
h.alt <- h.alt %>% select(COMID, contains("pn"))

#get COMIDs of occ datasets
h.alt.sub <- lapply(occ.sites, function(x) select(x, COMID, long, lat))
h.alt.sub <- bind_rows(h.alt.sub) %>% distinct()

h.alt.sub <- left_join(h.alt.sub, h.alt)

#fix NAs (no COMID in hydro alt dataset)
h.alt.sub <- h.alt.sub %>% mutate(index = row_number())

na.h.alt <- filter(h.alt.sub, is.na(pnSeasonal)) %>%
  st_as_sf(., coords = c("long", "lat"), crs = 4269, remove = FALSE) %>%
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
  
  na.h.alt$COMID_altNA <- as.character(nhd_near$COMID) #add COMID to data
  
  dist <- as.vector(st_distance(st_transform(na.h.alt, 5070), nhd_near, by_element = TRUE)) #get distance to nearest nhd strm + add to dataset
  na.h.alt$dist2strm_m_nhd_NA <- round(dist, 3)
  
  #rejoin w/ alt data and find if still NAs 
  tmp.na <- na.h.alt %>% 
    select(!contains(names(h.alt)), contains("COMID")) %>% 
    st_drop_geometry() %>%
    left_join(h.alt, by = c("COMID_altNA" = "COMID"))
  h.alt.sub <- bind_rows(filter(h.alt.sub, !index %in% tmp.na$index), tmp.na)
  
  #check for NAs still
  na.h.alt <- filter(h.alt.sub, is.na(pnSeasonal)) %>%
    st_as_sf(., coords = c("long", "lat"), crs = 4269, remove = FALSE) %>%
    st_transform(5070)
  na.coms <- c(na.coms, na.h.alt$COMID_altNA)
  
}

h.alt.sub <- h.alt.sub %>% 
  select(-index) %>%
  mutate(COMID_h_alt_source = case_when(is.na(COMID_altNA) ~ COMID, 
                                        T ~ COMID_altNA)) %>%
  select(-COMID_altNA) %>%
  relocate(contains("COMID"), contains("dist"))

write_csv(h.alt.sub, paste0(PATH, "/02_EnvDat/hydrologic_alteration_all_sites.csv"))

rm(h.alt, h.alt.sub, nhd.na, na.coms, na.h.alt, tmp.na, dist, nhd_near, nhd)

}

#Adjusted HIT metrics ----
if(get_adj_hit) {
  
#load in hit metrics + gage info
hit <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240130.csv"))
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv")) %>%
  select(-c(station_nm, site_tp_cd, coord_acy_cd, dec_coord_datum_cd, alt_va, alt_acy_va, alt_datum_cd, data_type_cd,
            parm_cd, stat_cd, ts_id, count_nu, period, int_high))

#intermittent streams and streams with zero flow days
adj.sites <- data.frame(site_no = hit$site_no, is_intermit = rep(NA, nrow(hit)), is_zero_flw = rep(NA, nrow(hit)), prop_zero_flw = rep(NA, nrow(hit)))

for(i in seq_len(nrow(adj.sites))) {
  site <- adj.sites$site_no[i]
  
  #check if flow type is intermittent
  adj.sites[i, ]$is_intermit <- g.info[g.info$site_no == site, ]$flw_type == "Int"
  
  #check for zero flow days
  #original raw flow dat
  f <- file.path(PATH, "02_EnvDat/raw_daily_flow", 
                 paste0(site, ifelse(file.exists(file.path(PATH, "02_EnvDat/raw_daily_flow", paste0(site, "_long.csv"))), #check if clipped file exists
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
  select(-all_of(adj_vars)) %>% #remove unadjusted version of variables
  rename_with(~gsub("_adj", "", .x))

hit <- bind_rows(
  hit[!hit$site_no %in% new.hit$site_no, names(hit)[!grepl("_adj", names(hit))]],
  new.hit
)

#save
write_csv(hit, paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_adjustments_applied.csv"))

rm(new.hit, adj.sites, adj_vars, g.info, hit)

}

# Env Change Over Time variables ----
## prep ----
### get date range of occurrence data ----
if(get_occ_date_range) {
 
  # load in dat with dates (if available)
  file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "inthigh_long", full.names = TRUE)
  occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
                                                           long = col_number(),
                                                           lat_new = col_number(),
                                                           long_new = col_number(),
                                                           taxa_count = col_number(),
                                                           .default = "c"))
  
  lapply(occ.list, function(x) {
    #filter to match wide data parameters
    x <- x %>%
      filter(dist2gage_m_15yr <= 15000 & dist2strm_m_flw <= 5000) %>% #make sure within 15km of a gage + 5 km classified flow line
      filter(!gage_no_15yr %in% c("06906800", "07020550", "07332500"))
    #combine all dates
    all.dates <- c(x$date_min, x$date_max, x$date)
    all.dates <- all.dates[!is.na(all.dates)]
    
    paste0(min(all.dates), " to ", max(all.dates))
    # return(all.dates)
  })
  
  # fish.dates <- all.dates[[2]]
  # fish.dates <- as.Date(fish.dates)
  #based on looking into the data a bit more, the ADEQ sites with a date of 1900 are likely wrong (looked at ADEQ website)
  
  #bugs: Aug. 1964 to Jan. 2022
  #fish: Jan. 1923 to Jun. 2021
   
}

##get bounding box for raster clipping ----
coords <- lapply(occ.sites, function(x) select(x, long, lat))
coords <- bind_rows(coords) %>% distinct()
coords <- st_as_sf(coords, coords = c("long", "lat"), crs = 4269)
coords <- st_bbox(coords)
#add some wiggle room
coords[c("xmax", "ymax")] <- coords[c("xmax", "ymax")] + 0.5
coords[c("xmin", "ymin")] <- coords[c("xmin", "ymin")] - 0.5

### get watersheds x nhd strms ----
huc <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/huc12_interior_highlands_crop.shp"))

#pair watersheds w/ comids
nhd <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/nhd_clip_state_eco.shp"))
nhd <- nhd[nhd$COMID %in% comids,]

if(file.exists(paste0(PATH, "/02_EnvDat/HUCs_NHDs/nhd_huc_intersection_info.csv"))) { #load intersection info if starting over
  
  huc_nhd_df <- read_csv(paste0(PATH, "/02_EnvDat/HUCs_NHDs/nhd_huc_intersection_info.csv"), col_types = cols(prop_in_huc = col_number(), .default = col_character()))
  
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

### functions ----
calc_var_diff <- function(x, #variable dataframe that has, at minimum, env variable + time
                          env_var = "air_temp", #column name for variable of interest
                          site_var = "huc12", #column name for individual sites to group by
                          var2grp = c("year", "season"), #group levels, right now this only works across years, not within a single year
                          norm_comp = TRUE, #calc deficit from normal, departure from normal, percent of normal (all averaged over time)
                          pct_change = TRUE, #average percent change across time for each lowest group
                          delta_var = TRUE, #difference in var from start to end of period, grouped by N years
                          num_yr_start_end = 5, #for difference from start to end of period, how many years to pull from
                          c_var_diff = TRUE, #cumulative change over length of time (essentially length of the trend line)
                          min_max_diff = TRUE, #difference between max value and min value across time period (magnitude ish)
                          growth_rate = TRUE #growth rate calculated from decomposed time series, not grouped by season or anything
) {
  # x <- air_temp #to test
  # add in site id column for grouping purposes
  var2grp <- c(site_var, var2grp)
  
  #set up dataframe for returning results
  all_res <- x[, names(x) %in% var2grp[!grepl("year", var2grp)]]
  all_res <- distinct(all_res)
  
  if(norm_comp) {
    #get normal = average value over all years ----
    normal <- x %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm = mean(.data[[env_var]], na.rm = TRUE))
    
    obs <- x %>%
      group_by(across(all_of(var2grp))) %>% 
      summarise(obs = mean(.data[[env_var]], na.rm = TRUE)) %>%
      ungroup() %>%
      left_join(., normal, by = names(.)[names(.) %in% names(normal)])
    
    #within years deficits ----
    def <- obs %>%
      mutate(def = (obs - norm) / norm) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_def = mean(def, na.rm = TRUE)) %>%
      ungroup()
    
    #departure from normal (obs - normal) ----
    depart <- obs %>%
      mutate(depart = obs - norm) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_depart = mean(depart, na.rm = TRUE)) %>%
      ungroup()
    
    #percentage of normal (obs/normal * 100) ----
    pct <- obs %>%
      mutate(pct = obs/norm*100) %>% 
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(norm_pct = mean(pct, na.rm = TRUE)) %>%
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
      summarise(mn_var = mean(.data[[env_var]])) %>%
      arrange(across(all_of(var2grp))) %>% 
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(diff_yr = year - lag(year),
             ann_pct_ch = ifelse(mn_var == 0 & lag(mn_var) == 0, 0,
                                 ifelse(lag(mn_var) == 0 & mn_var != 0, 100,
                                       ((mn_var / lag(mn_var))^(1 / diff_yr) - 1)*100)),
             ann_pct_ch = ifelse((mn_var - lag(mn_var)) >= 0, abs(ann_pct_ch), abs(ann_pct_ch)*-1)) %>% #second line to handle transitions from neg. to pos.
      # mutate(pct_ch = (mn_var / lag(mn_var) - 1) * 100,
      #        diff = ann_pct_ch - pct_ch) %>% 
      summarise(mn_pch = mean(ann_pct_ch, na.rm = TRUE))
    
    all_res <- all_res %>% left_join(., mn.pch, by = names(.)[names(.) %in% names(mn.pch)])
  }
  
  #diff between start and end of sample period ----
  if(delta_var) {
    dates <- unique(x$year)
    dates <- sort(dates)
    
    d_var <- x %>%
      mutate(time_point = case_when(year %in% head(dates, num_yr_start_end) ~ "first",
                                    year %in% tail(dates, num_yr_start_end) ~ "last")) %>%
      filter(!is.na(time_point)) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)])), time_point) %>%
      summarise(avg_var = mean(.data[[env_var]], na.rm = TRUE)) %>%
      arrange(across(all_of(var2grp[!grepl("year", var2grp)])), time_point) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(delta_var = avg_var - lag(avg_var)) %>%
      select(any_of(var2grp), delta_var) %>%
      filter(!is.na(delta_var))
    
    all_res <- all_res %>% left_join(., d_var,  by = names(.)[names(.) %in% names(d_var)])
  }
  
  
  #cumulative diff. of over time ----
  if(c_var_diff) {
    c_var <- x %>%
      group_by(across(all_of(var2grp))) %>%
      summarise(ssn_var = mean(.data[[env_var]])) %>% #one group value per year
      arrange(across(all_of(var2grp))) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      mutate(diff_yr = year - lag(year),
             cdelt_var = ssn_var - lag(ssn_var),
             adj_cdelt = abs(cdelt_var) / diff_yr) %>%
      summarise(cdelt_var = sum(abs(adj_cdelt), na.rm = TRUE))
    
    all_res <- all_res %>% left_join(., c_var,  by = names(.)[names(.) %in% names(c_var)])
  }
  
  #min/max diff ----
  if(min_max_diff) {
    mmd <- x %>%
      group_by(across(all_of(var2grp))) %>%
      summarise(grp_var = mean(.data[[env_var]])) %>%
      arrange(across(all_of(var2grp))) %>%
      group_by(across(all_of(var2grp[!grepl("year", var2grp)]))) %>%
      summarise(min = min(grp_var, na.rm = TRUE),
                max = max(grp_var, na.rm = TRUE),
                min_max_diff = max - min)
    
    all_res <- all_res %>% left_join(., mmd,  by = names(.)[names(.) %in% names(mmd)])
  }
  
  #growth rate from decomposed time series ----
  if(growth_rate) {
    ts.l <- split(x, x[[site_var]])
    
    ts.l <- lapply(ts.l, function(y) {
      y <- y %>% arrange(date)
      
      ts.gr <- ts(y[[env_var]], start = c(y[1, ]$year, y[1, ]$month), end = c(y[nrow(y), ]$year, y[nrow(y), ]$month), frequency = 12)
      
      m <- decompose(ts.gr, type = "multiplicative")
      # m <- decompose(ts.gr, type = "additive")
      trend <- na.omit(m$trend)
      trend <- diff(log(trend))
      
      #avg percent growth rate
      pct.avg.gr <- mean(trend) * 100
      
      res <- data.frame(site_var = unique(y[[site_var]]), pct_avg_gr = pct.avg.gr)
      
      names(res)[names(res) == "site_var"] <- site_var
      
      return(res)
    })
    
    ts.l <- bind_rows(ts.l)
    
    all_res <- all_res %>% left_join(., ts.l,  by = names(.)[names(.) %in% names(ts.l)])
  }
  
  names(all_res)[!names(all_res) %in% var2grp & !grepl("var", names(all_res))] <- paste0(names(all_res)[!names(all_res) %in% var2grp & !grepl("var", names(all_res))], "_", env_var)
  names(all_res)[grepl("var", names(all_res))] <- sub("var", env_var, names(all_res)[grepl("var", names(all_res))])
  
  return(all_res)
  
  #for testing:
  # rm(y, ts.gr, m, all_res, c_var, d_var, def, depart, mn.pch, normal, obs, pct, ts.l, x, c_var_diff, dates, delta_var, env_var, growth_rate, norm_comp, num_yr_start_end, pct_change, var2grp)
}

## air temp change ----
##note: prism data across almost all years were previously downloaded and cropped in 07_get_fstrm_temp_data.R
### summarize air temp within HUCs (Jan 1923 to Jan 2022) ----
if(get_huc_lvl_air_temp) {
  
  prism_files <- list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters"), full.names = TRUE)
  prism_files <- prism_files[grepl("tmean", prism_files)]
  
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
  
}

###taxa specific delta ----
if(summarize_air_temp) {
  
a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_huc_lvl.csv")) %>%
  #for seasonal trends:
  mutate(date = as.Date(paste(year, month, "01", sep = "-")),
         season = case_when(month %in% c("12", "01", "02") ~ "winter",
                            month %in% c("03", "04", "05") ~ "spring",
                            month %in% c("06", "07", "08") ~ "summer",
                            month %in% c("09", "10", "11") ~ "fall"))

####bugs: Aug. 1964 to Jan. 2022 ----
air_temp <- a_temp %>%
    filter(date >= as.Date("1964-08-01") & date <= as.Date("2022-01-31"))

bugs <- calc_var_diff(air_temp, env_var = "air_temp", site_var = "huc12",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

bugs <- bugs %>%
  mutate(across(where(is.numeric) & !matches("norm"), ~ round(.x, digits = 6)),
         across(matches("norm"), ~ round(.x, digits = 10)))

write_csv(bugs, paste0(PATH, "/02_EnvDat/air_temp_temporal_change_bugs.csv"))
  
###fish: Jan. 1923 to Jun. 2021 ----
air_temp <- a_temp %>%
  filter(date >= as.Date("1923-01-01") & date <= as.Date("2021-06-30"))

fish <- calc_var_diff(air_temp, env_var = "air_temp", site_var = "huc12",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

fish <- fish %>%
  mutate(across(where(is.numeric) & !matches("norm"), ~ round(.x, digits = 6)),
         across(matches("norm"), ~ round(.x, digits = 10)))

write_csv(fish, paste0(PATH, "/02_EnvDat/air_temp_temporal_change_fish.csv"))

### plots ----
# fish <- fish %>% arrange(huc12, season)
# bugs <- bugs %>% arrange(huc12, season)
# f.cor <- cor(bugs[, -c(1,2)])
# corrplot::corrplot(f.cor, method = "number")
# 
# ggplot() +
#   geom_point(data = fish, aes(delta_air_temp, cdelt_air_temp, color = season)) +
#   theme_minimal()
# 
# ggplot() +
#   geom_point(data = bugs, aes(delta_air_temp, cdelt_air_temp, color = season)) +
#   theme_minimal()
# 
# air_temp %>%
#   group_by(date, season) %>%
#   summarise(mn_temp = mean(air_temp)) %>%
#   ggplot() +
#   geom_line(aes(date, mn_temp, color = season)) +
#   theme_minimal()

#comparing diff change ratios
# air_temp %>%
#   #order is high neg, high pos, middle value
#   filter((huc12 == "071401020407" & season == "winter") | (huc12 == "102901070301" & season == "winter") | (huc12 == "111402010501" & season == "winter")) %>%
#   group_by(huc12) %>%
#   arrange(date) %>%
#   mutate(scaled_var = (air_temp - mean(air_temp)) / sd(air_temp)) %>%
#   ggplot(aes(date, scaled_var)) +
#   geom_line(aes(color = huc12), alpha = 0.3) +
#   geom_smooth(method = "loess", span = 0.7, aes(color = huc12)) + #change span to edit smoothing amount (higher is more)
#   theme_classic()

rm(a_temp, air_temp, bugs, f.cor, fish, output, date_min, date_max)

} ## end air temp section ##

## precipitation ----
### download data ----
if(download_prism_ppt) {
  
  library(prism)
  prism_set_dl_dir(paste0(PATH, "/02_EnvDat/raw_air_temp/prism"))
  
  get_prism_monthlys(type = "ppt",
                     years = 1923:2022,
                     mon = seq(1, 12),
                     keepZip = TRUE, #can remove the unzipped files after clipping
                     keep_pre81_months = TRUE)
  
}

### clip prism data to region ----
if(clip_prism_ppt) {
 
  prism.list <- list.dirs(paste0(PATH, "/02_EnvDat/raw_air_temp/prism"), full.names = FALSE, recursive = FALSE)
  prism.list <- prism.list[grepl("ppt", prism.list)]
  
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
  
  rm(prism.c, prism.rast, prism_files, prism.list, rast_name, i)
  
}

### avg monthly precip in occupied HUCs ----
if(get_huc_lvl_ppt) {

prism_files <- list.files(paste0(PATH, "/02_EnvDat/raw_air_temp/prism/clipped_rasters"), full.names = TRUE)
prism_files <- prism_files[grepl("ppt", prism_files)]

precip <- data.frame()
precip <- mclapply(prism_files, mc.cores = 8, function(x) { #~ 18 GB MEM at max using 8 cores
  date <- gsub(".*?(\\d{6}).*", "\\1", x)
  
  if(nchar(date) != 6) { return(NULL) }
  
  pr1 <- rast(x) #load in clipped raster
  
  #extract temperature values at all gage sites
  precip1 <- exact_extract(pr1, huc, fun = 'weighted_mean', weights = 'area', force_df = TRUE, progress = FALSE) #exactextractr::exact_extract is faster than terra::extract
  
  names(precip1) <- "precip"
  #add hucID and date to df
  precip1$year <- gsub("^(.{4}).*", "\\1", date)
  precip1$month <- gsub(".*(.{2})$", "\\1", date)
  precip1 <- bind_cols(st_drop_geometry(huc[, "huc12"]), precip1)
  
  return(precip1)
})

precip <- bind_rows(precip)
precip <- precip %>% arrange(huc12, year, month)

write_csv(precip, paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_precip_at_huc_lvl.csv"))

}

###taxa specific delta ----
if(summarize_ppt) {
  
precip <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_precip_at_huc_lvl.csv")) %>%
  #for seasonal trends:
  mutate(date = as.Date(paste(year, month, "01", sep = "-")),
         season = case_when(month %in% c("12", "01", "02") ~ "winter",
                            month %in% c("03", "04", "05") ~ "spring",
                            month %in% c("06", "07", "08") ~ "summer",
                            month %in% c("09", "10", "11") ~ "fall"))

####bugs: Aug. 1964 to Jan. 2022 ----
precip.sub <- precip %>%
  filter(date >= as.Date("1964-08-01") & date <= as.Date("2022-01-31"))

bugs <- calc_var_diff(precip.sub, env_var = "precip", site_var = "huc12",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

write_csv(bugs, paste0(PATH, "/02_EnvDat/precip_ssn_temporal_change_bugs.csv"))

bugs <- calc_var_diff(precip.sub, env_var = "precip", site_var = "huc12",
                      var2grp = c("year"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

write_csv(bugs, paste0(PATH, "/02_EnvDat/precip_all_temporal_change_bugs.csv"))

###fish: Jan. 1923 to Jun. 2021 ----
precip.sub <- precip %>%
  filter(date >= as.Date("1923-01-01") & date <= as.Date("2021-06-30"))

fish <- calc_var_diff(precip.sub, env_var = "precip", site_var = "huc12",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

write_csv(fish, paste0(PATH, "/02_EnvDat/precip_ssn_temporal_change_fish.csv"))

fish <- calc_var_diff(precip.sub, env_var = "precip", site_var = "huc12",
                      var2grp = c("year"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

write_csv(fish, paste0(PATH, "/02_EnvDat/precip_all_temporal_change_fish.csv"))

### plots ----
# ggplot(data = precip.sub %>% filter(huc12 == "111401090801" & season == "summer")) +
#   geom_line(aes(date, precip, color = season))
# 
# f.cor <- cor(fish[, -c(1,2)])
# corrplot::corrplot(f.cor, method = "number")
# 
# ggplot(data = bugs) +
#   # geom_point(aes(delta_precip, cdelt_precip, color = season)) +
#   geom_point(aes(norm_def_precip, mn_pch_precip)) +
#   theme_minimal()
# 
# ggplot(data = fish) +
#   # geom_point(aes(delta_precip, cdelt_precip, color = season)) +
#   # geom_point(aes(norm_def_precip, mn_pch_precip)) +
#   # geom_point(aes(norm_def_precip, norm_pct_precip)) +
#   geom_point(aes(norm_depart_precip, delta_precip)) +
#   theme_minimal()

rm(f.cor, fish, bugs, precip, precip.sub)

}

##stream temp ----
if(summarize_stream_temp) {
  
#stream temp dates are limited by the available dates of recorded water temp at usgs gages
pred_temp <- read_csv(paste0(PATH, "/02_EnvDat/predicted_stream_temps_monthly_envchange.csv"))
# gage_ids[!gage_ids %in% pred_temp$site_no] #double checking

pred_temp <- pred_temp[pred_temp$site_no %in% gage_ids, ]
pred_temp <- pred_temp %>%
  mutate(season = case_when(month %in% c("12", "1", "2") ~ "winter",
                            month %in% c("3", "4", "5") ~ "spring",
                            month %in% c("6", "7", "8") ~ "summer",
                            month %in% c("9", "10", "11") ~ "fall"))

### bugs: Aug. 1964 to Jan. 2022 ----
water_temp <- pred_temp %>%
  filter(date >= as.Date("1964-08-01") & date <= as.Date("2022-01-31"))

bugs <- calc_var_diff(water_temp, env_var = "water_temp", site_var = "site_no",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

bugs <- bugs %>%
  mutate(across(where(is.numeric) & !matches("norm"), ~ round(.x, digits = 6)))

write_csv(bugs, paste0(PATH, "/02_EnvDat/stream_temp_temporal_change_bugs.csv"))

# water_temp %>%
#   group_by(date, season) %>%
#   summarise(mn_temp = mean(water_temp)) %>%
#   ggplot() +
#   geom_line(aes(date, mn_temp, color = season)) +
#   theme_minimal()

## fish: Jan. 1923 to Jun. 2021 ----
water_temp <- pred_temp %>%
  filter(date >= as.Date("1923-01-01") & date <= as.Date("2021-06-30"))

fish <- calc_var_diff(water_temp, env_var = "water_temp", site_var = "site_no",
                      var2grp = c("year", "season"),
                      norm_comp = TRUE, pct_change = TRUE, 
                      delta_var = TRUE, num_yr_start_end = 5,
                      c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = TRUE)

fish <- fish %>%
  mutate(across(where(is.numeric) & !matches("norm"), ~ round(.x, digits = 6)))

write_csv(fish, paste0(PATH, "/02_EnvDat/stream_temp_temporal_change_fish.csv"))

## plots ----
# f.cor <- cor(bugs[, -c(1,2)])
# corrplot::corrplot(f.cor, method = "number")
# 
# water_temp %>%
#   filter(site_no == "06923500") %>%
#   group_by(year, season) %>%
#   summarise(mn_temp = mean(water_temp)) %>%
#   ggplot() +
#   geom_line(aes(year, mn_temp, color = season)) +
#   scale_y_continuous(limits = c(0, 30))
#   
# ggplot(data = fish) +
#   # geom_point(aes(delta_water_temp, cdelt_water_temp, color = season)) +
#   # geom_point(aes(delta_water_temp, pct_avg_gr_water_temp, color = season)) +
#   geom_point(aes(norm_def_water_temp, cdelt_water_temp, color = season)) +
#   theme_minimal()
# 
# ggplot(data = bugs) +
#   # geom_point(aes(delta_water_temp, cdelt_water_temp, color = season)) +
#   # geom_point(aes(delta_water_temp, pct_avg_gr_water_temp, color = season)) +
#   geom_point(aes(delta_water_temp, pct_avg_gr_water_temp, color = season)) +
#   theme_minimal()

rm(bugs, fish, pred_temp, water_temp)

}

##land cover change ----
###download nlcd data ----
if(download_nlcd) {
  
library(FedData)

#make download template
clip_box <- matrix(c(coords[["xmin"]], coords[["ymin"]],
                     coords[["xmax"]], coords[["ymin"]],
                     coords[["xmax"]], coords[["ymax"]],
                     coords[["xmin"]], coords[["ymax"]],
                     coords[["xmin"]], coords[["ymin"]]), 
                     ncol = 2, byrow = TRUE)
clip_box <- st_polygon(list(clip_box)) %>%
  st_sfc()
clip_box <- st_sf(clip_box, crs = 4269)

#break down box into downloadable sizes
multi_clip_box <- st_make_grid(clip_box, n = c(5, 5))

#set downloadable years
download_years <- c(2001, 2004, 2006, 2008, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)
#dataset types
data_types <- c("landcover", "impervious", "canopy")

#for loop to download section x year
for(i in seq_along(multi_clip_box)) {
  
  for(n_year in download_years) {
    
    for(dat.t in data_types) {
      
      if(dat.t %in% c("landcover", "impervious") & n_year %in% c(2012, 2014, 2015, 2017, 2018, 2020)) { next } #not available in these years
      if(dat.t == "canopy" & n_year %in% c(2001, 2004, 2006, 2008)) { next } #not available in these years
      
      if(dat.t == "landcover") { data_name <- "Land_Cover" }
      if(dat.t == "impervious") { data_name <- "Impervious" }
      if(dat.t == "canopy") { data_name <- "Tree_Canopy" }

      if(file.exists(paste0(PATH, "/02_EnvDat/NLCD/box", i, "_NLCD_", data_name, "_", n_year, ".tif"))) { next }
      
      tmp <- get_nlcd(
        template = multi_clip_box[i],
        dataset = dat.t,
        year = n_year,
        label = paste0("box", i),
        extraction.dir = file.path(PATH, "02_EnvDat", "NLCD")
      )
      
      message(dat.t, ": ", n_year, " for box", i, " downloaded")
      rm(data_name)
      Sys.sleep(0.5)
      
    }
    
  }
  
  message(round(i/length(multi_clip_box)*100, digits = 2), "% complete")
  
}

rm(clip_box, multi_clip_box, tmp, download_years, i, n_year)

} ## end nlcd download section ##

###summarize lc by huc ----
if(get_huc_lvl_nlcd) {
#match nlcd crs
huc.tr <- st_transform(huc, 5070)

#imperviousness and tree canopy
dat_years <- c(2001, 2004, 2006, 2008, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021)
data_types <- c("Impervious", "Tree_Canopy")

dat_summ <- expand.grid(huc12 = huc_ids) # set up dataframe to get results
for(data_name in data_types) { #for loop for mem
  
  message("Starting ", data_name, " summarizing")
  
  for(n_year in dat_years) {
    
    if(data_name == "Impervious" & n_year %in% c(2012, 2014, 2015, 2017, 2018, 2020)) { next } #not available in these years
    if(data_name == "Tree_Canopy" & n_year %in% c(2001, 2004, 2006, 2008)) { next } #not available in these years
    
    file_list <- list.files(file.path(PATH, "02_EnvDat", "NLCD"), pattern = paste0(data_name, "_", n_year, ".tif"), full.names = TRUE)
    file_list <- file_list[!grepl(".xml", nlcd_list)]
    
    tmp <- vrt(file_list) #load the tiled rasters as one layer
    
    #extract temperature values at all gage sites
    dat1 <- exact_extract(tmp, huc.tr, function(df) {
      df %>%
        mutate(value_frac = value * coverage_fraction) %>%
        group_by(huc12) %>%
        summarise(pct_cov = sum(value_frac, na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE))
    }, summarize_df = TRUE, include_cols = 'huc12', progress = FALSE)
    
    dat1 <- dat1 %>%
      mutate(pct_cov = round(pct_cov, digits = 6))
    
    names(dat1)[names(dat1) == "pct_cov"] <- paste(data_name, "pct", n_year, sep = "_")
    
    dat_summ <- dat_summ %>%
      left_join(., dat1, by = "huc12")
    
    message(round(which(n_year == dat_years)/length(dat_years)*100, digits = 2), "% complete.")
    
  }
  
}

write_csv(dat_summ, paste0(PATH, "/02_EnvDat/NLCD/imp_treecanopy_yearly_at_huc_lvl.csv"))

#land cover
nlcd_years <-  c(2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019, 2021)

nlcd_summ <- expand.grid(huc12 = huc_ids, nlcd_class = nlcd_colors()$ID) # set up dataframe to get results
for(n_year in nlcd_years) {
  nlcd_list <- list.files(file.path(PATH, "02_EnvDat", "NLCD"), pattern = paste0("Land_Cover_", n_year, ".tif"), full.names = TRUE)
  nlcd_list <- nlcd_list[!grepl(".xml", nlcd_list)]
  
  tmp <- vrt(nlcd_list) #load in list of rasters for specific date, join together
  
  #extract temperature values at all gage sites
  nlcd1 <- exact_extract(tmp, huc.tr, function(df) {
    df %>%
      mutate(frac_total = coverage_fraction / sum(coverage_fraction)) %>%
      group_by(huc12, value) %>%
      summarize(freq = sum(frac_total))
  }, summarize_df = TRUE, include_cols = 'huc12', progress = FALSE)
  
  nlcd1 <- nlcd1 %>%
    mutate(freq = round(freq*100, digits = 6))
  
  names(nlcd1)[names(nlcd1) == "freq"] <- paste0("pct_cov_", n_year)
  names(nlcd1)[names(nlcd1) == "value"] <- "nlcd_class"
  
  nlcd_summ <- nlcd_summ %>%
    left_join(., nlcd1, by = names(.)[names(.) %in% names(nlcd1)])
  
  message(round(which(n_year == nlcd_years)/length(nlcd_years)*100, digits = 2), "% complete.")
}

nlcd_summ <- nlcd_summ %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.x), 0, .x)))

write_csv(nlcd_summ, paste0(PATH, "/02_EnvDat/NLCD/nlcd_yearly_at_huc_lvl.csv"))

rm(nlcd_list, tmp, nlcd1, nlcd_summ, dat1, huc.tr, dat.t, dat_years, data_types, nlcd_years, n_year, file_list)
}

###summarize change over time ----
if(summarize_nlcd) {
#no taxonomic specific dates bc of temporal extents
nlcd_summ <- read_csv(paste0(PATH, "/02_EnvDat/NLCD/nlcd_yearly_at_huc_lvl.csv"))
dat_summ <- read_csv(paste0(PATH, "/02_EnvDat/NLCD/imp_treecanopy_yearly_at_huc_lvl.csv"))

tmp1 <- nlcd_summ %>%
  pivot_longer(contains("pct_cov"), names_to = "year", values_to = "nlcd_pct") %>%
  mutate(year = gsub("pct_cov_", "", year),
         year = as.numeric(year),
         nlcd_class = as.factor(nlcd_class))

nlcd_delt <- calc_var_diff(tmp1, env_var = "nlcd_pct", site_var = "huc12",
                           var2grp = c("year", "nlcd_class"),
                           norm_comp = TRUE, pct_change = TRUE, 
                           delta_var = TRUE, num_yr_start_end = 2,
                           c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = FALSE)

write_csv(nlcd_delt, paste0(PATH, "/02_EnvDat/nlcd_temporal_change_alltax.csv"))

#prep canopy cover and imperviousness
tmp2 <- dat_summ %>%
  pivot_longer(contains("pct_"), names_to = "type_year", values_to = "dat_pct") %>%
  mutate(year = gsub(".*_pct_", "", type_year),
         year = as.numeric(year),
         data_type = gsub("_pct_.*$", "", type_year)) %>%
  select(-type_year)

#canopy cover
tcc <- tmp2 %>% 
  filter(data_type == "Tree_Canopy") %>%
  rename(tcc_pct = dat_pct) %>%
  select(-data_type)

tcc_delt <- calc_var_diff(tcc, env_var = "tcc_pct", site_var = "huc12",
                          var2grp = c("year"),
                          norm_comp = TRUE, pct_change = TRUE, 
                          delta_var = TRUE, num_yr_start_end = 2,
                          c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = FALSE)

write_csv(tcc_delt, paste0(PATH, "/02_EnvDat/tcc_temporal_change_alltax.csv"))

#imp
imp <- tmp2 %>% 
  filter(data_type == "Impervious") %>%
  rename(imp_pct = dat_pct) %>%
  select(-data_type)

imp_delt <- calc_var_diff(imp, env_var = "imp_pct", site_var = "huc12",
                          var2grp = c("year"),
                          norm_comp = TRUE, pct_change = TRUE, 
                          delta_var = TRUE, num_yr_start_end = 2,
                          c_var_diff = TRUE, min_max_diff = TRUE, growth_rate = FALSE)

write_csv(imp_delt, paste0(PATH, "/02_EnvDat/impervious_temporal_change_alltax.csv"))


####plots ----
# f.cor <- cor(nlcd_delt[, -c(1,2)])
# corrplot::corrplot(f.cor, method = "number")
# 
# imp %>%
#   filter(huc12 == "111101030303") %>%
#   # group_by(year) %>%
#   # summarise(imp_pct = median(imp_pct)) %>%
#   ggplot() +
#   geom_line(aes(year, imp_pct)) +
#   theme(legend.position = "none")
# tcc %>%
#   filter(huc12 == "110702090501") %>%
#   # group_by(year) %>%
#   # summarise(tcc_pct = median(tcc_pct)) %>%
#   ggplot() +
#   geom_line(aes(year, tcc_pct)) +
#   theme(legend.position = "none")

# plot(huc.tr[huc.tr$huc12 == "111401070407",]$geometry)
# plot(tmp, add = TRUE)
# nlcd_summ %>%
#   group_by(nlcd_class) %>%
#   summarise(across(where(is.numeric), mean)) %>%
#   pivot_longer(contains("pct_cov"), names_to = "year", values_to = "pct_cov") %>%
#   mutate(year = gsub("pct_cov_", "", year),
#          year = as.numeric(year),
#          nlcd_class = as.factor(nlcd_class)) %>%
#   ggplot() +
#   geom_line(aes(year, pct_cov, color = nlcd_class))

rm(nlcd_delt, tmp1, nlcd_summ, tmp2, imp, imp_delt, tcc, tcc_delt, dat_summ)
}

# Variable Summ. + Plots ----
# source("Scripts/XX_colors.R")
# 
# flow.link <- lapply(occ.sites, function(x) {
#   x %>%
#     select(COMID, gage_no_15yr, flw_type) %>%
#     rename(site_no = gage_no_15yr)
# })
# flow.link <- bind_rows(flow.link)
# flow.link <- distinct(flow.link)
# 
# flow.link <- left_join(flow.link, select(huc_nhd_df, -comid_idx), by = c("COMID" = "comid"))
# 
# ## stream cat ----
# stream_cat <- read_csv(paste0(PATH, "/02_EnvDat/StreamCat/stream_cat_vars_all_sites.csv"), col_types = cols(COMID = col_character(), .default = col_number()))
# stream_cat <- lapply(occ.sites, function(x) {
#   x %>%
#     select(COMID, flw_type) %>%
#     left_join(., stream_cat)
# })
# 
# stream_cat[[1]]$taxa <- "bug"
# stream_cat[[2]]$taxa <- "fish"
# 
# stream_cat <- bind_rows(stream_cat)
# 
# stream_cat %>%
#   select(-taxa) %>%
#   distinct() %>%
#   ggplot() +
#   geom_density(aes(CLAYWS, fill = flw_type), alpha = 0.6) +
#   scale_y_continuous(expand = c(0,0)) +
#   # geom_boxplot(aes(flw_type, TMAX8110WS, fill = flw_type)) +
#   scale_fill_manual(values = flow.pal) +
#   theme_classic()
# 
# ## hydro alt ----
# hydro_alt <- read_csv(paste0(PATH, "/02_EnvDat/hydrologic_alteration_all_sites.csv"))
# 
# ## HIT ----
# hit <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_adjustments_applied.csv"))
# # hit.orig <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_metrics_20240130.csv"))
# 
# hit <- lapply(occ.sites, function(x) {
#   x %>%
#     select(gage_no_15yr, flw_type) %>%
#     rename(site_no = gage_no_15yr) %>%
#     left_join(., hit)
# })
# 
# hit[[1]]$taxa <- "bug"
# hit[[2]]$taxa <- "fish"
# 
# hit <- bind_rows(hit)
# 
# hit %>%
#   # select(-taxa) %>%
#   filter(taxa == "fish") %>%
#   distinct() %>%
#   ggplot() +
#   geom_density(aes(dh1, fill = flw_type), alpha = 0.6) +
#   scale_y_continuous(expand = c(0,0)) +
#   # geom_boxplot(aes(flw_type, TMAX8110WS, fill = flw_type)) +
#   scale_fill_manual(values = flow.pal) +
#   theme_classic()
