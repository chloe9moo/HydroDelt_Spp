##PREDICT STRM TEMP

# modified from the methods of Middaugh et al. 2016 (Ecology of Freshwater Fish)

#STEPS:
#1. Calculate least-squares linear regression for each site to predict stream temperature from air temperature
#2. Predict stream temperature at GF sites using matching flow or matching county regressions (??)
#3. Summarize stream temps for each site

library(tidyverse); library(sf)

PATH <- getwd()

#1. Calculate least-squares linear regression for each site to predict stream temperature from air temperature ----
a_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_air_temp/prism_monthly_temp_at_strm_gage_sites.csv"), 
                   col_types = cols(site_no = col_character(), Lat = col_character(), Long = col_character(),
                                    .default = col_number()))

w_temp <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_temp_20240314.csv"),
                   col_types = cols(site_no = col_character(), Lat = col_character(), Long = col_character(),
                   .default = col_number())) %>%
  rename(water_temp = mn_monthly_temp) %>%
  filter(ct >= 20 & !is.na(Lat)) #remove sites with fewer than 20 days in the monthly mean and with missing coord info

c_temp <- left_join(w_temp, a_temp)
# write_csv(c_temp, paste0(PATH, "/02_EnvDat/combined_water_air_temp.csv"))

c_temp <- split(c_temp, c_temp$site_no)
# lapply(c_temp, function(site){
#   cat(unique(site$site_no), "has", sum(is.na(site$water_temp))/nrow(site) * 100, "percent NAs in water temp\n")
# })

r_list <- lapply(c_temp, function(df) {
  
  if(all(is.na(df$water_temp))) { message("No data for df", unique(df$site_no), "\n"); return(NULL) }
  
  df <- df %>% arrange(year, month) %>% filter(!is.na(water_temp))
  lsr <- lm(water_temp ~ air_temp, data = df)
  lsr$site_no <- unique(df$site_no)
  
  return(lsr)
})
saveRDS(r_list, paste0(PATH, "/02_EnvDat/fine_scale_temp_regressions.rds"))

##plot to look at regression ----
plot_temp_regression <- function(model_list = r_list, site_no, comb_df_list = c_temp) {
  if(!is.null(model_list)) {
    p1 <- ggplot(data = model_list[[site_no]]$model, aes(air_temp, water_temp)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, lwd = 0.7) +
      scale_color_viridis_c() +
      theme_classic() +
      ggtitle(paste0("site ", site_no))
  } else { p1 <- NULL }
  
  if(!is.null(comb_df_list)) {
    new_df <- comb_df_list[[site_no]] %>%
      pivot_longer(cols = contains("_temp"), names_to = "type", values_to = "temp") %>%
      arrange(year, month, type) %>%
      mutate(date = as.Date(paste(year, month, "01", sep = "-")))

    p2 <- ggplot(data = new_df, aes(x = date, y = temp, group = type, color = type)) +
      geom_path() +
      # scale_x_date(date_breaks = "6 month", date_labels = "%Y-%m") +
      theme_classic()
  } else { p2 <-  NULL }
  
  ggpubr::ggarrange(p1, p2, nrow = 2)
  
}

plot_temp_regression(r_list, "06818000", c_temp)

ggsave(paste0(PATH, "/99_figures/ex_fn_strm_temp_regression.png"), dpi = 300, height = 5, width = 6, bg = "white")
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

rm(l1, i)

# predict stream temp at ungaged sites ----
g.info <- read_csv(paste0(PATH, "/02_EnvDat/raw_stream_temp/comb_water_monthly_location_info.csv")) %>%
  filter(!is.na(Lat)) %>%
  st_as_sf(., coords = c("Long", "Lat"), remove = FALSE, crs = 4269)

#find which models are useable 
#looking for model performance but also models built on at least a year of data
pred_list <- r_list[names(r_list) %in% lm_res[lm_res$sig < 0.05 & lm_res$R2 > 0.8 & !is.na(lm_res$sig), ]$site]
keep <- sapply(pred_list, function(x) nrow(x$model) >= 12) #1 yr of data
pred_list <- pred_list[keep]

#find which sites need predicted stream temp, get available strm dates
hit.info <- read_csv(paste0(PATH, "/02_EnvDat/usgs_gage_hit_inthigh_allinfo.csv")) %>%
  st_as_sf(., coords = c("dec_long_va", "dec_lat_va"), remove = FALSE, crs = 4269)

#get distance to nearest water temp gage within same flow type
lm.sites <- g.info[g.info$site_no %in% names(pred_list),]
nr.site <- st_nearest_feature(hit.info, lm.sites, check_crs = TRUE) #find nearest model source
lm_near <- lm.sites[nr.site,]
hit.info$strm_lm_site_no <- lm_near$site_no #add site no to match sites

##get distance
dist <- as.vector(st_distance(hit.info, lm_near, by_element = TRUE)) #get distance to nearest nhd strm
hit.info$dist2strm_temp_m <- round(dist, digits = 2)

write_csv(hit.info, paste0(PATH, "/02_EnvDat/raw_stream_temp/usgs_gage_info_with_strm_temp_info.csv"))

rm(lm_near, nr.site, dist, keep)

##predict stream temp using nearest gage
pred_temp <- data.frame()
for(i in seq_len(nrow(hit.info))) {
  
  site_no <- hit.info[i, ]$site_no #gage of interest
  lm_site_no <- hit.info[i, ]$strm_lm_site_no #matched gage
  
  lm <- r_list[names(r_list) == lm_site_no][[1]] #regression from nearest gage
  
  lm.w_temp <- c_temp[[lm_site_no]] %>% #get dates of water temp from regression
    mutate(date = as.Date(paste(year, month, "1", sep = "-")))
  
  sub.a_temp <- a_temp[a_temp$site_no == site_no, ] %>%
    mutate(date = as.Date(paste(year, month, "1", sep = "-"))) %>%
    filter(date <= max(lm.w_temp$date) & date >= min(lm.w_temp$date)) %>%
    arrange(date)
  
  sub.a_temp$water_temp <- predict(lm, newdata = sub.a_temp["air_temp"]) #predict stream temp
  
  sub.a_temp$lm_site_no <- lm_site_no #save site model is based on (JIC)
  
  sub.a_temp$temp_type <- "predicted" #for comparisons later (JIC)
  
  sub.a_temp <- relocate(sub.a_temp, site_no, date, Lat, Long, water_temp, air_temp)
  
  if(site_no == lm_site_no) { #modify data if observed water temp is available
    
    sub.a_temp <- sub.a_temp %>% 
      left_join(., lm.w_temp %>% 
                  select(date, water_temp) %>% 
                  rename(orig_water_temp = water_temp))%>%
      mutate(water_temp = case_when(is.na(orig_water_temp) ~ water_temp,
                                    T ~ orig_water_temp),
             temp_type = "observed") %>%
      select(-orig_water_temp)
    
  }
  
  pred_temp <- bind_rows(pred_temp, sub.a_temp)
  
}
rm(site_no, lm_site_no, lm, sub.a_temp, lm.w_temp)

write_csv(pred_temp, paste0(PATH, "/02_EnvDat/predicted_stream_temps_monthly.csv"))

# summarize stream temp values ----
#occ data goes from 1900 to 2022
pred_temp <- read_csv(paste0(PATH, "/02_EnvDat/predicted_stream_temps_monthly.csv"))

tmp <- pred_temp %>%
  #prep df
  mutate(date = as.Date(date),
         season = case_when(month %in% c(12, 1, 2) ~ "winter",
                            month %in% c(3, 4, 5) ~ "spring",
                            month %in% c(6, 7, 8) ~ "summer",
                            month %in% c(9, 10, 11) ~ "fall")) %>%
  #calc summ stats for each month
  group_by(site_no, month) %>%
  mutate(across(water_temp, list(
    mn = ~mean(.x, na.rm = T),
    med = ~median(.x, na.rm = T),
    max = ~max(.x, na.rm = T),
    min = ~min(.x, na.rm = T)
  ), .names = "{.fn}_mnth_wtemp")) %>%
  ungroup() %>%
  #calc summ stats for each season
  group_by(site_no, season) %>%
  mutate(across(water_temp, list(
    mn = ~mean(.x, na.rm = T),
    med = ~median(.x, na.rm = T),
    max = ~max(.x, na.rm = T),
    min = ~min(.x, na.rm = T),
    cv = ~(sd(.x, na.rm = T) / mean(.x, na.rm = T))), 
    .names = "{.fn}_ssn_wtemp")) %>%
  ungroup() %>%
  #get avg annual CV
  group_by(site_no, year) %>%
  mutate(cv_yr_wtemp = sd(water_temp, na.rm = T) / mean(water_temp, na.rm = T)) %>%
  ungroup()

mnth <- tmp %>%
  select(site_no, month, contains("mnth")) %>%
  distinct() %>%
  pivot_wider(id_cols = site_no, 
              names_from = month,
              values_from = contains("mnth"),
              names_glue = "{.value}_{month}")
ssn <- tmp %>%
  select(site_no, season, contains("ssn")) %>%
  distinct() %>%
  pivot_wider(id_cols = site_no, 
              names_from = season,
              values_from = contains("ssn"),
              names_glue = "{.value}_{season}") 
cv <- tmp %>%
  select(site_no, contains("yr")) %>%
  distinct() %>%
  group_by(site_no) %>%
  summarize(ann_wtemp_cv = mean(cv_yr_wtemp, na.rm = TRUE))

summ.temp <- left_join(mnth, ssn) %>% left_join(., cv)

#thermal sensitivity from slopes
tmp <- split(pred_temp, pred_temp$site_no)

slope.l <- lapply(tmp, function(x) {
  lm <- lm(water_temp ~ air_temp, data = x)
  slope <- coef(lm)[2]
  return(slope)
})

slope.l <- data.frame(site_no =  names(slope.l), therm_sens = unlist(slope.l))

summ.temp <- summ.temp %>%
  left_join(., slope.l)

write_csv(summ.temp, paste0(PATH, "/02_EnvDat/predicted_stream_temps_summ_all.csv"))
