##CALCULATE FINE STREAM TEMP FROM REGRESSIONS

# following the methods of Middaugh et al. 2016 (Ecology of Freshwater Fish)

#STEPS:
#1. Calculate least-squares linear regression for each site to predict stream temperature from air temperature
#2. Predict stream temperature at GF sites using matching flow or matching county regressions (??)

library(tidyverse); library(sf)

PATH <- getwd()

#1. Calculate least-squares linear regression for each site to predict stream temperature from air temperature ----
a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/noaa_air_monthly_temp_countylevel_waterIDS.csv")) %>% rename(air_temp = mn_monthly_temp) %>% select(-mn_max_monthly_temp, -mn_min_monthly_temp)
w_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_temp_20231025.csv")) %>%
  rename(water_temp = mn_monthly_temp) %>%
  filter(ct >= 20) %>%
  mutate(date = with(., sprintf("%d-%02d", year, month)),
         p_date = as.Date(paste0(date, "-01")))
c_temp <- left_join(a_temp, w_temp, by = c("site_no", "p_date", "date", "year", "month")) %>%
  relocate(site_no, fips, date, air_temp, water_temp)
# write_csv(c_temp, paste0(PATH, "/02_EnvDat/combined_air_water_temp.csv"))

c_temp <- split(c_temp, c_temp$site_no)
# lapply(c_temp, function(site){
#   cat(unique(site$site_no), "has", sum(is.na(site$water_temp))/nrow(site) * 100, "percent NAs in water temp\n")
# })

r_list <- lapply(c_temp, function(site) {
  
  if(all(is.na(site$water_temp))) { cat("No data for site", unique(site$site_no), "\n"); return(NULL) }
  
  site <- site %>% arrange(date) %>% filter(!is.na(water_temp))
  lsr <- lm(water_temp ~ air_temp, data = site)
  lsr$site_no <- unique(site$site_no)
  
  return(lsr)
})
saveRDS(r_list, paste0(PATH, "/02_EnvDat/fine_scale_temp_regressions.rds"))

##plot to look at regression ----
ggplot(data = r_list$`07150500`$model, aes(air_temp, water_temp)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, lwd = 0.7) +
  scale_color_viridis_c() +
  theme_classic() +
  ggtitle("")
ggsave(paste0(PATH, "/99_figures/ex_fn_strm_temp_regression.png"), dpi = 300, height = 4, width = 4, bg = "white")
# plot(r_list$`07198000`$model$water_temp, r_list$`07198000`$residuals)

##1a. Make table of linear regression results to examine ----
lm_res <- data.frame()
for(i in 1:length(r_list)) {
  
  if(is.null(r_list[[i]])) {
    
    l1 <- data.frame(site = names(r_list)[[i]], intercept = NA, coeff = NA, R2 = NA, flow_type = NA)
  
  } else {
    
    l1 <- data.frame(
      site = r_list[[i]]$site_no,
      intercept = coef(r_list[[i]])[[1]],
      coeff = coef(r_list[[i]])[[2]],
      R2 = summary(r_list[[i]])$r.squared,
      flow_type = NA)
    
    if(summary(r_list[[i]])$df[2] == 0) { #for the models that didn't work (low count)
      l1$sig <- NA
    } else {
      l1$sig <- summary(r_list[[i]])$coefficients["air_temp", "Pr(>|t|)"]
    }
  }
  
  lm_res <- bind_rows(lm_res, l1)
}

lm_res <- lm_res %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3)))
write_csv(lm_res, paste0(PATH, "/02_EnvDat/fine_scale_temp_regression_results.csv"))

rm(l1, i, a_temp)

#2. predict stream temp ----
r_list <- readRDS(paste0(PATH, "/02_EnvDat/fine_scale_temp_regressions.rds")) #regressions
lm_res <- read_csv(paste0(PATH, "/02_EnvDat/fine_scale_temp_regression_results.csv")) #regression results
g.xy <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_prevHITsites_countyinfo.csv")) #county info on sites (usgs gage and other dat)
water_sites <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/all_water_temp_locations_info.csv")) %>% #location info on strm temp sites
  st_as_sf(., coords = c("Long", "Lat"), crs = 4269)
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_all_prevHITsites_info.csv")) %>% #target gage location info
  select(site_no, dec_long_va, dec_lat_va) %>%
  distinct() %>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269)
us.air.temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/usgs_air_monthly_temp_countylevel.csv")) #compiled air temp from climate viewer
hit <- read_csv(paste0(PATH, "/02_EnvDat/Updated_HIT_4.23.a.csv")) #hit data (orig. target gages)

#which counties are in the study area?
hlnd <- st_as_sf(maps::map("county", c("arkansas", "oklahoma", "missouri"), fill = T, plot = F)) %>% #EPSG 4269 = NAD83
  st_transform(., crs = 4269)
eco <- read_sf(paste0(PATH, "/02_EnvDat/study_extent_shp/ecoreg_l3_interior_highlands_crop.shp")) %>%
  st_union()
sf_use_s2(FALSE)
ct.in.eco <- hlnd %>% st_filter(., eco)
ct.in.eco <- bind_rows(ct.in.eco, hlnd %>% filter(ID %in% c("arkansas,desha","missouri,st louis city")))

#filter out NULL models or poor performing models
r_list <- r_list[!sapply(r_list, is.null)] #remove NULLs
r_list <- r_list[names(r_list) %in% lm_res[lm_res$sig < 0.05 & lm_res$R2 > 0.8, ]$site]

pred_temps <- data.frame()
for (i in 1:nrow(hit)) {
  
  id <- hit[i, ]$STAID_0
  fips <- g.xy[grepl(id, g.xy$site_no), ]$fips
  
  if(!any(g.xy[g.xy$fips == fips, ]$county %in% ct.in.eco$ID)) { #skip sites outside the study extent, not worth calculating right now
    cat("County not in study extent, skipping to next site...\n"); next
  }
  
  n.at <- us.air.temp %>% filter(fips_code == fips) #pull out usgs climate viewer data for that county
  
  ##opt 2.1: if the site has a specific model run for it ----
  if(any(names(r_list) == id)) { 
    
    lm <- r_list[names(r_list) == id][[1]] #pull out model
    model_type <- "Orig. model for target site"
    
  } else {
     
  #search for other models in that county with a model
    sub.id <- g.xy[grepl(fips, g.xy$fips), ]$site_no
    sub.id <- sub.id[sub.id != id]
    
  #opt 2.2: #if there is a model alternative found in the same county ----
    if(any(names(r_list) %in% sub.id)) {
      if(length(sub.id) == 1) { #if there is one model option in the same county, pull it out
        
        lm <- r_list[names(r_list) == sub.id][[1]]
        model_type <- "Same county, only one model alternative"
        
      } else { #otherwise ...
        
        lm <- r_list[names(r_list) %in% sub.id]
        
        if(length(lm) == 1) { # it was actually only one alternative
          
          lm <- lm[[1]]
          model_type <- "Same county, only one model alternative"
          
        } else {
          
          #get all xy coords of interest
          target.station <- g.info %>% filter(site_no == id)
          surr.stations <- water_sites %>% filter(site_no %in% names(lm))
          
          #find closest station
          dist <- st_distance(target.station, surr.stations)
          y <- which.min(dist)
          y <- names(lm)[[y]]
          
          lm <- r_list[names(r_list) == y][[1]]
          model_type <- "Same county, selected nearest site from multiple"
          
        }
      } 
  #opt 2.3: if no model alternative is found in that county ----
    } else { 
  
      ##opt 2.3.1: search for surrounding county models ----
    t.ct <- hlnd %>%
      filter(ID %in% g.xy[g.xy$fips == fips, ]$county) #get target county shape
    
    surr.ct <- hlnd[st_touches(hlnd, t.ct, sparse = FALSE),] #get surrounding counties
    
    near.sites <- g.xy[g.xy$county %in% surr.ct$ID, ]$site_no #get sites from surrounding counties
    
    if(any(names(r_list) %in% near.sites)) {
      
      lm <- r_list[names(r_list) %in% near.sites]
      
      if(length(lm) == 1) { #check to see if more than one surrounding county is pulled out
        
        lm <- r_list[names(r_list) %in% near.sites][[1]]
        model_type <- "Surrounding county model, only one model alternative"
        
      } else {
        
        #get all xy coords of interest
        target.station <- g.info %>% filter(site_no == id)
        surr.stations <- water_sites %>% filter(site_no %in% names(lm))
        
        #find closest station
        dist <- st_distance(target.station, surr.stations)
        y <- which.min(dist)
        y <- names(lm)[[y]]
        
        lm <- r_list[names(r_list) == y][[1]]
        model_type <- "Surrounding county model, selected nearest site from multiple"
        
      }
      
    } else {
      
      #get all xy coords of interest
      target.station <- g.info %>% filter(site_no == id)
      surr.stations <- water_sites %>% filter(site_no %in% names(r_list))
      
      #find closest station
      dist <- st_distance(target.station, surr.stations)
      y <- which.min(dist)
      y <- surr.stations[y, ]$site_no
      
      lm <- r_list[names(r_list) == y][[1]]
      model_type <- "Surrounding county model, selected nearest site from multiple"  
      
      }
    }
  } #end site doesn't have specific model for it

  #predict for both 4.5 and 8.5 climate models
  pred.dat <- n.at %>% select(mn_temp_8.5) %>% rename(air_temp = mn_temp_8.5)
  n.at$pred_strm_temp_8.5 <- predict(lm, newdata = pred.dat)
  
  pred.dat <- n.at %>% select(mn_temp_4.5) %>% rename(air_temp = mn_temp_4.5)
  n.at$pred_strm_temp_4.5 <- predict(lm, newdata = pred.dat)
  
  n.at <- n.at %>% mutate(across(contains("pred_strm"), ~ round(.x, digits = 3))) #round to same digits as climate data
  n.at$site_no <- id
  
  n.at <- n.at %>% 
    select(date, fips_code, site_no, contains("pred_strm")) %>%
    mutate(model_from = model_type)
  
  pred_temps <- bind_rows(pred_temps, n.at)
}

write_csv(pred_temps, paste0(PATH, "/02_EnvDat/predicted_stream_temps_monthly.csv"))

#plotting for examples ----
water_sites <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/all_water_temp_locations_info.csv")) %>%
  st_as_sf(., coords = c("Long", "Lat"), crs = 4269)
g.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_all_prevHITsites_info.csv")) %>%
  select(site_no, dec_long_va, dec_lat_va) %>%
  distinct() %>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), crs = 4269)

ct.highlight <- hlnd %>% filter(., ID %in% g.xy[g.xy$fips == fips, ]$county)

target.station <- g.info %>% filter(site_no == id)
near.stations <- water_sites %>% 
  filter(site_no %in% y) #%>%
  mutate(type = ifelse(site_no == "06906800", "No strm temp data to model from", "Has strm temp data"))

ggplot() +
  geom_sf(data = hlnd) +
  geom_sf(data = ct.in.eco, fill = "#CC8C97") +
  # geom_sf(data = surr.ct, fill = "black") +
  geom_sf(data = ct.highlight, fill = "red") +
  geom_sf(data = surr.stations, color = "white") +
  geom_sf(data = near.stations, color = "blue")+
  geom_sf(data = target.station, color = "black") +
  scale_color_manual(values = c("black", "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(paste0(PATH, "/99_figures/ex_fn_strm_temp_noclosemodeloptions.png"), dpi = 300, width = 8, height = 5, bg = "white")
