## OBTAIN + PREP ENV VARIABLES FOR GF MODELS ##

library(tidyverse); library(sf); library(terra)
options(readr.show_col_types = FALSE)

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
#get date range of occurrence data
#load in dat with dates (if available)
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
# }) 

#bugs: Aug. 1964 to Jan. 2022
#fish: Jan. 1900 to Jun. 2021

YOU ARE HERE FOR MONDAY!!!!

#from NOAA - https://water.weather.gov/precip/ (some ideas here)
#precip deficit (???)
#departure from normal? (not sure how this would work if we need 1 value for each site)
#average departure from normal? could also do amount of variation??


#stream temp change

#air temp change

#precip change



#land cover change

##ADD HDI VARIABLES FROM FOX & MAGOULICK 2022?







t <- nhd %>%
  mutate(missing_c = ifelse(COMID %in% m.c, "wrong", "right")) %>%
  group_by(missing_c) %>%
  summarise(across(where(is.numeric), mean))




