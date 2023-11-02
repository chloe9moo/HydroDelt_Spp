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
write_csv(c_temp, paste0(PATH, "/02_EnvDat/combined_air_water_temp.csv"))

c_temp <- split(c_temp, c_temp$site_no)
# lapply(c_temp, function(site){
#   cat(unique(site$site_no), "has", sum(is.na(site$water_temp))/nrow(site) * 100, "percent NAs in water temp\n")
# })

r_list <- lapply(c_temp, function(site) {
  
  if(all(is.na(site$water_temp))) { cat("No data for site\n", unique(site$site_no)); return(NULL) }
  
  site <- site %>% arrange(date) %>% filter(!is.na(water_temp))
  lsr <- lm(site$water_temp ~ site$air_temp)
  lsr$site_no <- unique(site$site_no)
  
  return(lsr)
})
saveRDS(r_list, paste0(PATH, "/02_EnvDat/fine_scale_temp_regressions.rds"))

# ggplot(data = r_list$`310`$model, aes(`site$air_temp`, `site$water_temp`)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   scale_color_viridis_c() +
#   theme_classic()
# plot(r_list$`07198000`$model$`site$water_temp`, r_list$`07198000`$residuals)

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
      l1$sig <- summary(r_list[[i]])$coefficients["site$air_temp", "Pr(>|t|)"]
    }
  }
  
  lm_res <- bind_rows(lm_res, l1)
}

lm_res <- lm_res %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3)))
write_csv(lm_res, paste0(PATH, "/02_EnvDat/fine_scale_temp_regression_results.csv"))

rm(l1, i, a_temp)






