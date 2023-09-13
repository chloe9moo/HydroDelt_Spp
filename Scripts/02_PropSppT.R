## CALC. PROP. SPECIES W/ THRESHOLD
# by C.E.Moore V.
# 
# OBJ: Calculate proportion of species with a threshold

library(tidyverse); library(gradientForest)
# library(parallel)

PATH <- getwd()
# source(paste0(PATH, "/Scripts/01a_FindThreshFunc.R"))

file.name <- "full_thresh_gf.bugs.23.GW_full_cat.csv" #<<which threshold file to analyze // single model runs
# file.list <- list.files(paste0(PATH, "/11_Thresholds"), pattern = "full_thresh", full.names = F) #parallel / multi model runs

#<< what proportion of overall cum imp curve must the threshold be to count (C2) + proportion of presences (C3, C4)? unit = percent
t_prop_val <- 25 

# n.cores <- detectCores()/2
# mclapply(file.list, mc.cores = n.cores, function(file.name) { #start parallel run, ~6 min for 18 models, up to 14.5 GB RAM

t <- read_csv(paste0(PATH, "/11_Thresholds/", file.name)) 

gf.mod <- readRDS(paste0(PATH, "/10_GFOutput/", gsub(".csv", "", gsub("full_thresh_", "", file.name)), ".rds")) #<<original gf model for associated threshold results
spp_imp <- as.data.frame(importance(gf.mod, type = "Species")) %>% rownames_to_column("species") #pull out species performance
colnames(spp_imp)[2] <- "rel_err"
full.in <- cbind(gf.mod$Y, gf.mod$X) #pull out full predictor + response data and combine (for later)

iVars <- unique(t$env_var)
nspecies <- length(unique(t$species))

p.t.df <- data.frame(env_var = iVars, C_1 = NA)

#Criteria 1. Spp with ANY threshold ----
for(var in iVars) {
  tmp <- t %>% 
    filter(env_var == var) %>%
    group_by(species) %>%
    filter(sum(slope, na.rm = T) != 0) #get rid of spp definitely no threshold (no change)
  p.t.df[p.t.df$env_var == var, "C_1"] <- (length(unique(tmp$species)) / nspecies) * 100
  rm(tmp)
}

cat("\nCriteria 1 calculated..\n")

#Criteria 2. Threshold = T based on proportion of overall cum imp change ----
tmp <- t %>%
  select(species, env_var, t_max_prop, t_m2_prop) %>%
  distinct() %>%
  mutate(n_t = if_else(t_m2_prop > t_max_prop & !is.na(t_m2_prop), t_m2_prop, t_max_prop)) %>% #take larger change threshold
  group_by(env_var) %>%
  summarise(C_2 = sum(n_t > t_prop_val, na.rm = T)) %>%
  mutate(C_2 = (C_2 / nspecies) * 100)
p.t.df <- left_join(p.t.df, tmp)

cat("\nCriteria 2 calculated..\n")

#Criteria 3. calculate proportion of presence absence on either side of threshold ----
spp.list <- as.list(unique(t$species))
tmp <- t %>%
  select(species, env_var, t_max_mn, t_max_prop, t_m2_mn, t_m2_prop) %>%
  distinct() %>%
  mutate(t_mn = if_else(t_m2_prop > t_max_prop & !is.na(t_m2_prop), t_m2_mn, t_max_mn))

c3.df <- data.frame()
for(var in iVars) {
  
  pa_distr <- lapply(spp.list, function(spp) {
    # cat(spp, "\n") #for debug
    tmp2 <- tmp %>% filter(species == spp & env_var == var)
    
    if(is.na(tmp2$t_mn)) {
      pa <- data.frame(species = spp, P_A = factor(c(0,1)), ct_var_greater_t=NA, ct_var_less_t=NA, p_diff = NA, Q1 = NA, pa_diff = NA, Q2 = NA)
      return(pa)
    }
    
    pa <- full.in %>% 
      select(one_of(spp, var)) %>%
      rename(P_A = all_of(spp)) %>% #fix column name for later.
      group_by(P_A) %>%
      mutate(ct = n()) %>%
      reframe(ct_var_greater_t = (sum(.data[[var]] > unique(tmp2$t_mn), na.rm = T) / ct)*100,
              ct_var_less_t = (sum(.data[[var]] < unique(tmp2$t_mn), na.rm = T) / ct)*100) %>%
      distinct() %>%
      arrange(P_A) %>%
#1. Is the threshold capturing majority of presences on one side?
      mutate(p_diff = abs(ct_var_greater_t - ct_var_less_t),
             Q1 = if_else(p_diff > t_prop_val, TRUE, FALSE)) %>%
#2. Is the threshold capturing a larger portion of presences than absences on that side?
      mutate(max = if_else(p_diff == 0, 0, pmax(ct_var_greater_t, ct_var_less_t))) %>% #which side is greater?
      mutate(max = case_when(max == ct_var_greater_t ~ "ct_var_greater_t",
                             max == ct_var_less_t ~ "ct_var_less_t",
                             max == 0 ~ "tie")) %>%
      pivot_longer(contains("ct_"))
    
    col <- unique( pa[pa$P_A==1, ]$max)
    if(col != "tie") {
      tmp2 <- pa %>% filter(name == col)
      Q2 <- tmp2[tmp2$P_A == 1, ]$value > tmp2[tmp2$P_A == 0, ]$value
      pa_diff <- abs(tmp2[tmp2$P_A == 1, ]$value - tmp2[tmp2$P_A == 0, ]$value)
    } else {
      Q2 <- TRUE
      pa_diff <- 50 - pa[which(50 > pa[pa$P_A == 0, ]$value), ]$value
      if(length(pa_diff) == 0) { Q2 <- FALSE; pa_diff <- 0 }
    }
      
    pa <- pa %>% pivot_wider() %>%
      mutate(species = spp, pa_diff = pa_diff, Q2 = Q2) %>% select(-max) %>%
      relocate(species, P_A, ct_var_greater_t, ct_var_less_t, p_diff, Q1, pa_diff, Q2)
  
    return(pa)
  })
  
  pa_distr <- bind_rows(pa_distr)
  pa_distr$env_var <- var
  
  c3.df <- rbind(c3.df, pa_distr)
  cat(var, "criteria 3 complete..\n")
}

write_csv(c3.df, paste0(PATH, "/11_Thresholds/P_A_distrib_", gsub("full_thresh_", "", file.name)))

p.t.df <- c3.df %>%
  filter(!is.na(ct_var_greater_t) & P_A == 1) %>%
  filter(pa_diff > t_prop_val) %>%
  group_by(env_var) %>%
  summarise(C_3 = sum(Q1 & Q2)) %>%
  mutate(C_3 = (C_3 / nspecies) * 100) %>%
  right_join(., p.t.df) %>%
  relocate(env_var, C_1, C_2, C_3)

cat("\nCriteria 3 calculated..\n")

#Criteria 4. Criteria 2 and Criteria 3 together ----
## threshold is greater than 25% of overall cumulative importance change AND can distinguish between presence / absence
tmp <- t %>%
  select(species, env_var, t_max_prop, t_m2_prop) %>%
  distinct() %>%
  mutate(n_t = if_else(t_m2_prop > t_max_prop & !is.na(t_m2_prop), t_m2_prop, t_max_prop)) %>% #take larger change threshold
  filter(n_t > t_prop_val & !is.na(n_t))

tmp2 <- c3.df %>%
  filter(!is.na(ct_var_greater_t) & P_A == 1) %>%
  filter(pa_diff > t_prop_val & Q1 == T & Q2 == T)

for(var in iVars) {
  n <- length(intersect(tmp[tmp$env_var == var, ]$species, #species meeting criteria 2
                        tmp2[tmp2$env_var == var, ]$species)) #species meeting criteria 3
  p.t.df[p.t.df$env_var == var, "C_4"] <- (n / nspecies) * 100
}

cat("\nCriteria 4 calculated..\n")

#save prop csv ----
p.t.df <- p.t.df %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
write_csv(p.t.df, paste0(PATH, "/11_Thresholds/PropThresh_", gsub("full_thresh_", "", file.name)))

# }) #end parallel

# Criteria X. Spp above ratio dens split importance:dens obs ----
#from Chen et al. 2023 -- think this isn't relevant, it's regarding density not cumulative importance (ie splits density plots, not cum imp plots)
#make table of ratios
# r.df <- data.frame(env_var = iVars, ratio = NA)
# for (var in iVars) {
#   r.df[r.df$env_var == var, ]$ratio <- get_dens_ratio(gf.mod, var)
# }
#find species with threshold density ratio > 1

  
