### GRADIENT FORESTS ###

library(tidyverse); library(gradientForest)
options(readr.show_col_types = FALSE)

PATH <- getwd()

#read in data ----
#bio 
# file.list <- list.files(paste0(PATH, "/01_BioDat"), pattern = "_wide_", full.names = TRUE)
# occ.list <- lapply(file.list, read_csv, col_types = cols(lat = col_number(),
#                                                          long = col_number(),
#                                                          site_id = col_character(),
#                                                          COMID = col_character(),
#                                                          gage_no_15yr = col_character(),
#                                                          dist2gage_m_15yr = col_number(),
#                                                          dist2strm_m_flw = col_number()))

#trait
f <- list.files(paste0(PATH, "/20_Traits"), pattern = "site_x_trait", full.names = TRUE)
site.list <- lapply(f, function(x) {
  read_csv(x) %>%
    rename_with(~ gsub("-", "_", .x))
})

#env dat
env.dat <- read_csv(paste0(PATH, "/02_EnvDat/env_dat_combined_allsites.csv")) %>%
  rename(gage_lat = dec_lat_va, gage_long = dec_long_va) %>%
  select(-huc_cd, -contains("date"), -contains("dist"), -hit_ttl_yrs, -flw_name,
         -contains("mnth"), -contains("4.5"), -contains("ssn_max"), -contains("ssn_min"), -contains("ssn_med"))

#get variable list ----
vars.all <- names(env.dat)[!names(env.dat) %in% c("site_no", "COMID", "flw_type", "gage_class")]

vars.all <- data.frame(var = c(vars.all, "lat", "long"))

#add in group vars for separate runs
vars.all <- vars.all %>%
  mutate(fine_cat = case_when(grepl("pn|HDI", var) ~ "hydro_alt",
                              grepl("mnth|ssn|ann", var) ~ "stream_temp", 
                              grepl("lat|long|drain", var) ~ "spatial",
                              grepl("amplitude|ar|dh|dl|fh|fl|lam|ma|mh|ml|phase|ra|ta|th|tl", var) ~ "HIT"),
         course_cat = case_when(grepl("hydro_alt|stream_temp|HIT", fine_cat) ~ "hydrology",
                                grepl("spatial", fine_cat) ~ "spatial"), 
         run = case_when(var %in% vars[[1]] ~ TRUE, 
                         var %in% c("lat", "long", "pnSeasonal") ~ TRUE,
                         var %in% paste0("pn", str_to_upper(vars[[1]])) ~ TRUE,
                         T ~ FALSE)) #THIS IS FROM 10_GRADIENT_FORESTS_ORIGDAT.R!!!!!!!

#get gage types within each flow type
table(env.dat$gage_class, env.dat$flw_type)

#run GF models ----
#make function for repeated runs
gf_sub <- function(full.data, #biological data + env data dataframe
                   species.list, #species of interest (in case less than the full dataframe)
                   env.variable.list = vars.all, #variables of interest; previously pulled out HIT, LULC, and full run
                   flow.type = NA, #flow type to run, set as NA if all
                   gage.type = NA, #gage type to run, set as NA if all
                   gf.classification = TRUE, #classification or regression? likely don't change ever for our purposes
                   #gradientForest options:
                   ntree = 999,
                   transform = NULL,
                   compact = TRUE,
                   trace = TRUE,  
                   nbin = 201,
                   corr.threshold = 0.5) {
  
  ##remove sites with < 5 spp ----
  full.data <- full.data %>% 
    mutate(loc_tot = rowSums(select(., any_of(species.list)))) %>%
    filter(loc_tot > 5) %>%
    select(-loc_tot)
  
  ##subset by flow + reference ----
  #flow
  if(!any(grepl(flow.type, c("Int", "RO", "GW"), ignore.case = TRUE) | is.na(flow.type))) { stop("flow type not recognized.") }
  if(!is.na(flow.type)) { full.data <- full.data[grepl(flow.type, full.data$flw_type, ignore.case = TRUE), ] }
  
  #reference
  if(!any(grepl(gage.type, c("ref", "non-ref"), ignore.case = TRUE) | is.na(gage.type))) { stop("gage type not recognized.") }
  if(!is.na(gage.type)) { full.data <- full.data[grepl(paste0("^", gage.type, "$"), full.data$type, ignore.case = TRUE), ] }
  
  ##pull out env vars ----
  ##- remove NA env var rows
  full.data <- full.data %>% filter(!if_any(all_of(env.variable.list), is.na))
  env_col <- full.data[, grepl(paste(paste0("^", env.variable.list, "$"), collapse = "|"), colnames(full.data))]
  if(any(!env.variable.list %in% colnames(env_col))) { 
    warning(paste("The following environmental variables were not in the dataframe:", paste(env.variable.list[!env.variable.list %in% colnames(env_col)], collapse = ", "))) 
  }
  ## remove constant env vars
  col_rem <- sapply(env_col, function(col) all(col == col[1]))
  env_col <- env_col[, !col_rem]
  
  ##pull out species ----
  spp_col <- full.data[, grepl(paste(species.list, collapse = "|"), colnames(full.data))]
  ##- remove species w/ < 10 records
  spp_col <- spp_col[, colSums(spp_col) >= 10]
  # ##- remove species w/ < 5% of max collection records?? ##not sure they did this in the end...
  # spp_col <- spp_col[, colSums(spp_col) >= max(colSums(spp_col))*0.05]
  ##- remove species present at every site
  spp_col <- spp_col[, colSums(spp_col) != nrow(spp_col)]
  
  ##set model parameters ----
  n.sites <- nrow(spp_col)
  n.spp <- ncol(spp_col)
  l <- floor(log2(n.sites * 0.368/2))
  
  ##run as classification or regression? ----
  if(gf.classification) {
    spp_col <- spp_col %>% mutate(across(where(is.numeric), factor))
  }
  
  ##run model ----
  message("Running GF on ", n.spp, " species and ", n.sites, " sites\nUsing ", flow.type, " flow and ", gage.type, " gages..\n")
  set.seed(31)
  
  gf.out <- gradientForest(cbind(env_col, spp_col),
                           predictor.vars = colnames(env_col), response.vars = colnames(spp_col),
                           ntree = ntree, transform = transform, compact = compact, trace = trace,  
                           nbin = nbin, maxLevel = l, corr.threshold = corr.threshold)
  
  return(gf.out)
}

#testing ----
dir.create(paste0(PATH, "/10_GFOutput/", gsub("-", "_", Sys.Date())))
#bug tests
comb.dat <- left_join(site.list[[1]], env.dat %>% select(-flw_type, -COMID), by = c("gage_no_15yr" = "site_no"))

test <- gf_sub(full.data = comb.dat,
               species.list = names(site.list[[1]] %>% select(-c(lat, long, COMID, contains("flw"), site_id, contains("gage")))),
               env.variable.list = vars.all[vars.all$run == TRUE, ]$var)
saveRDS(test, paste0(PATH, "/10_GFOutput/2024_02_21/gf_bug_traits_allflw_allgage_all.rds"))

#break out flows
flow.opt <- c("Int", "RO", "GW")

for(i in seq_along(site.list)) {
  
  if(i == 2) { bio <- "fish"} else { bio <- "bugs" } #set taxa
  
  comb.dat <- left_join(site.list[[i]], env.dat %>% select(-flw_type, -COMID), by = c("gage_no_15yr" = "site_no"))
  spp <- names(site.list[[i]] %>% select(-c(lat, long, COMID, contains("flw"), site_id, contains("gage"))))
  
  for(f in flow.opt) {
    
    cat("\nRunning", bio, ":", f, "flow...\n")
   
    file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_traits_", f, "_all.rds")
    if(file.exists(file.name)) { cat("Model exists, next.\n"); next }
    
    gf.out <- NULL
    
    gf.out <- gf_sub(full.data = comb.dat, #biological data + env data dataframe
                     species.list = spp, #species of interest (in case less than the full dataframe)
                     env.variable.list = vars.all[vars.all$run == TRUE, ]$var, #variables of interest; previously pulled out HIT, LULC, and full run
                     flow.type = f, #flow type to run, set as NA if all
                     gage.type = NA #gage type to run, set as NA if all
    )
    
    saveRDS(gf.out, file = file.name) #saving bc it would be massive to store all of them in memory
    
    Sys.sleep(1)
    
  }
}


#run several groups of models at once: ----
#models to run: 
#fish + bug (2)
# all[[1]], all[[2]]
#all flows, int, ro, gw
flow.opt <- c("Int", "RO", "GW")
#all gages, ref, non-ref
gage.opt <- c(NA, "ref", "non-ref")
#all vars, HIT, LULC, H_ALT
var.opt <- list(vars.all$var, 
                vars.all[vars.all$var_type == "HIT", ]$var, 
                vars.all[vars.all$var_type == "LULC", ]$var, 
                vars.all[vars.all$var_type == "H_ALT", ]$var)

#run models
dir.create(paste0(PATH, "/10_GFOutput/", gsub("-", "_", Sys.Date()))) #date directory for future comparisons

for(i in 1:length(all)) {
  for(f in flow.opt) {
    for(g in gage.opt) {
      for(v in 1:length(var.opt)) {
        
        if(v == 1) { var <- "all" }
        if(v == 2) { var <- "HIT" }
        if(v == 3) { var <- "LULC" }
        if(v == 4) { var <- "HyALT" }
        if(i == 1) { bio <- "fish"} else { bio <- "bugs" }
        if(is.na(g)) { ga <- "all" } else { ga <- g }
        
        cat("\nRunning", bio, ":", f, "flow,", ga, "gage,", var, "variable type...\n")
        # file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_", f, "_", ga, "gage_", var, ".rds")
        file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_", f, "_", ga, "gage_", var, ".rds") #forgot x y vars
        if(file.exists(file.name)) { cat("Model exists, next.\n"); next }
        
        gf.out <- NULL

        gf.out <- gf_sub(full.data = all[[i]], #biological data + env data dataframe
                         species.list = spp_list$col_new, #species of interest (in case less than the full dataframe)
                         env.variable.list = var.opt[[v]], #variables of interest; previously pulled out HIT, LULC, and full run
                         flow.type = f, #flow type to run, set as NA if all
                         gage.type = g #gage type to run, set as NA if all
        )

        saveRDS(gf.out, file = file.name) #saving bc it would be massive to store all of them in memory

        Sys.sleep(1)

      }
    }
  }
}

