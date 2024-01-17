## FLOW CALC + ESTIMATION ##

library(tidyverse); library(sf)
library(dataRetrieval)
options(readr.show_col_types = FALSE)

PATH <- getwd()

# pull daily flow ----
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
eco <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "ecoreg_l3_interior_highlands_crop") #int. highlands boundary
# huc <- st_read(paste0(PATH, "/02_EnvDat/study_extent_shp/"), layer = "huc10_interior_highlands_crop")

flow.gage <- flow.gage %>%
  mutate(period = as.Date(end_date) - as.Date(begin_date)) %>%
  #15 years of data
  filter(period >= 15 * 365) %>%
  #make spatial for int high filter
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269) #NAD83
  
ih.flow.gage <- flow.gage %>%
  st_filter(., st_buffer(eco, dist = 2000)) #use int. highland eco region with 2 km buffer (for wiggle room)

#add flag to full ARMOOK set and save (for later filtering as necessary)
flow.gage <- flow.gage %>%
  mutate(int_high = ifelse(site_no %in% ih.flow.gage$site_no, TRUE, FALSE))
write_csv(st_drop_geometry(flow.gage), paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"))

#check
# ggplot() +
#   geom_sf(data = flow.gage, aes(color = int_high)) +
# #   geom_sf(data = ih.flow.gage, color = "blue") +
#   geom_sf(data = eco, fill = NA, color = "red") +
# #   # geom_sf(data = huc, fill = NA, color = "blue") +
#   theme_minimal()

rm(tmp, eco)

## retrieve flow data for filtered gages ----
# ih.flow.gage <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_flow_ARMOOK_info.csv"), 
#                          col_types = cols(site_no = col_character())) %>% #make sure sites load in correctly
#   filter(int_high == TRUE)

site_nm <- unique(ih.flow.gage$site_no)

for(i in 1:length(site_nm)) { #pull flow data and save
  nm <- site_nm[[i]]
  
  if(file.exists(paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))) { #in case it fails at any point
    message("Flow dat for ", nm, " already exists."); next
  } else {
    
    tmp.dat <- readNWISdv(siteNumbers = nm,
                          parameterCd = "00060",
                          statCd = "00003",
                          startDate = min(ih.flow.gage[ih.flow.gage$site_no == nm, ]$begin_date),
                          endDate = max(ih.flow.gage[ih.flow.gage$site_no == nm, ]$end_date))
    
    tmp.dat <- tmp.dat %>%
      mutate(day_diff = Date - lag(Date)) #for removing sites with > 7 
    message("Site ", nm, " has ", max(tmp.dat$day_diff, na.rm = TRUE), " day(s) at most between records.")
    
    write_csv(tmp.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", nm, ".csv"))
    message("Flow dat for site ", nm, " downloaded.")
    
  }
}

#figure out if there is a 15 yr period with < 7 contiguous missing days at each site, remove those without
site_check <- data.frame(site_no = site_nm, action = NA, max_diff = NA, years = NA)

for(i in 1:nrow(site_check)) {
  tmp.dat <- read_csv(paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check[i, ]$site_no, ".csv"), col_types = cols(site_no = col_character()))
  
  yrs <- as.numeric(tmp.dat[nrow(tmp.dat), ]$Date - tmp.dat[1,]$Date) / 365 #possible # of years
  
  if(max(tmp.dat$day_diff, na.rm = TRUE) <= 7) { #if whole dataset has < 7 missing consecutive days, keep it all/no change
    site_check[i, ]$action <- "keep all"
    site_check[i, ]$max_diff <- max(tmp.dat$day_diff, na.rm = TRUE)
    site_check[i, ]$years <- yrs
    message("Site ", site_check[i, ]$site_no, " meets all time requirements, no change to data.")
    
  } else {
    
    while(yrs > 15) { #clip until less than 15 potential years to pull from
      
      ind1 <- which.max(tmp.dat$day_diff) #when is the first big time break
      
      #find the next > 7 day break
      ind2 <- which(tmp.dat$day_diff[(ind1 + 1):nrow(tmp.dat)] >= 7)[1] + ind1 
      
      if(is.na(ind2)) { ind2 <- nrow(tmp.dat) } #if NA, then go to the end of the df (no other time breaks)
      
      yrs <- as.numeric(tmp.dat[ind2, ]$Date - tmp.dat[ind1, ]$Date) / 365 #number of years in new section
      
      if(yrs > 15) { #if the section is greater than 15 years, clip to period and save

        tmp.dat <- tmp.dat[ind1:ind2, ]
        tmp.dat[1, ]$day_diff <- NA
        
        write_csv(tmp.dat, paste0(PATH, "/02_EnvDat/raw_daily_flow/", site_check[i, ]$site_no, "_v2.csv"))
        message("Site ", site_check[i, ]$site_no, " had > 7 missing contig. days, clipped and update saved.")
        
        site_check[i, ]$action <- "clip"
        site_check[i, ]$max_diff <- max(tmp.dat$day_diff, na.rm = TRUE)
        site_check[i, ]$years <- yrs
        
        break #break while loop, found data meeting criteria
        
      } else {
        
        #remove that period and try again
        tmp.dat <- tmp.dat[-c(ind1:ind2), ] %>% 
          mutate(day_diff = Date - lag(Date))
        
        yrs <- as.numeric(tmp.dat[nrow(tmp.dat), ]$Date - tmp.dat[1,]$Date) / 365 #new possible # of years
        
      } #end if else - clipping
    } #end while loop
    
    if(yrs < 15) {
      message("Site ", site_check[i, ]$site_no, " too many missing consecutive days, not enough years.")
      break
    }
    
  }
}



# 1. Find USGS gages with daily stream flow measurments in the region
#     1a. probably pull all gages in ARMOOK, filter to region, then pull actual flow data
#     1b. Needs:
#       15 complete years of data
#       < 7 days of contiguous missing data
#       range standardize hydrologic data for gages (by dividing each metric by their absolute values to transform metrics to comparable scales?)
# 
# 2. Use EflowStats to calc HIT metrics (171 metrics) + 6

