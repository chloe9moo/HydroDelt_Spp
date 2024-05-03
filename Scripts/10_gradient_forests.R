### GRADIENT FORESTS ###

library(tidyverse); library(gradientForest)
options(readr.show_col_types = FALSE)

PATH <- getwd()

#read + prep data ----
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
fl <- list.files(paste0(PATH, "/20_Traits"), pattern = "site_x_trait_.*_abundance", full.names = TRUE)
site.list <- lapply(fl, function(x) {
  x <- read_csv(x) %>%
    select(-contains("dist2")) %>%
    rename(site_no = gage_no_15yr)
  names(x) <- gsub("-", "_", names(x))
  return(x)
})

species.list <- names(site.list[[3]])[!names(site.list[[3]]) %in% c("lat", "long", "COMID", "flw_type", "site_no", "dist2gage_m_15yr", "dist2strm_m_flw", "site_id")]

#env dat
env.info <- read_csv(paste0(PATH, "/02_EnvDat/environmental_variable_info.csv"))
env.info <- env.info[env.info$include_baseline == "yes" & !is.na(env.info$include_baseline), ] #only include previously selected variables

var.names <- env.info$variable_col_name

env.list <- lapply(unique(env.info$variable_source)[unique(env.info$variable_source) != "taxa dataset"], function(x) {
  file.name <- paste0(PATH, x)
  env.file <- read_csv(file.name)
  env.file <- env.file[, names(env.file) %in% c(var.names, "COMID", "site_no")]
  return(env.file)
})

#combine bio + env dat
for(i in seq_along(env.list)) {
  site.list[[3]] <- left_join(site.list[[3]], env.list[[i]])
}

#get variable list ----
# vars.all <- names(env.dat)[!names(env.dat) %in% c("site_no", "COMID", "flw_type", "gage_class")]
# 
# vars.all <- data.frame(var = c(vars.all, "lat", "long"))
# 
# #add in group vars for separate runs
# vars.all <- vars.all %>%
#   mutate(fine_cat = case_when(grepl("pn|HDI", var) ~ "hydro_alt",
#                               grepl("mnth|ssn|ann", var) ~ "stream_temp", 
#                               grepl("lat|long|drain", var) ~ "spatial",
#                               grepl("amplitude|ar|dh|dl|fh|fl|lam|ma|mh|ml|phase|ra|ta|th|tl", var) ~ "HIT"),
#          course_cat = case_when(grepl("hydro_alt|stream_temp|HIT", fine_cat) ~ "hydrology",
#                                 grepl("spatial", fine_cat) ~ "spatial"), 
#          run = case_when(var %in% vars[[1]] ~ TRUE, 
#                          var %in% c("lat", "long", "pnSeasonal") ~ TRUE,
#                          var %in% paste0("pn", str_to_upper(vars[[1]])) ~ TRUE,
#                          T ~ FALSE)) #THIS IS FROM 10_GRADIENT_FORESTS_ORIGDAT.R!!!!!!!

#get gage types within each flow type
# table(env.dat$gage_class, env.dat$flw_type)

#run GF models ----
#make function for repeated runs
gf_sub <- function(bio.env.dat = site.list[[3]], #biological data + env data dataframe
                   bio.type = "trait", #if biological vs. trait data, either "trait" or "taxonomic"
                   species.names = species.list, #species of interest (in case less than the full dataframe)
                   env.variable.names = var.names, #variables of interest; previously pulled out HIT, LULC, and full run
                   flow.type = NA, #flow type to run, set as NA if all
                   gage.type = NA, #gage type to run, set as NA if all
                   gf.classification = FALSE, #classification or regression? when to use each?
                   #gradientForest options:
                   ntree = 999,
                   transform = NULL,
                   compact = TRUE,
                   trace = TRUE,  
                   nbin = 201,
                   corr.threshold = 0.5) {
  
  ##remove sites with < 5 spp ----
  #WILL NEED TO FIGURE OUT A WAY TO ACCOUNT FOR THIS IN TRAIT DATA - PROB. NEED TO REMOVE BEFORE OR SOMETHING
  if(bio.type == "taxonomic") {
    bio.env.dat <-   bio.env.dat %>% 
      mutate(loc_tot = rowSums(select(., any_of(species.list)))) %>%
      filter(loc_tot > 5) %>%
      select(-loc_tot)
  }
  
  ##subset by flow + reference ----
  #flow
  if(!any(grepl(flow.type, c("Int", "RO", "GW"), ignore.case = TRUE) | is.na(flow.type))) { stop("flow type not recognized.") }
  if(!is.na(flow.type)) { bio.env.dat <- bio.env.dat[grepl(flow.type, bio.env.dat$flw_type, ignore.case = TRUE), ] }
  
  #reference
  if(!any(grepl(gage.type, c("ref", "non-ref"), ignore.case = TRUE) | is.na(gage.type))) { stop("gage type not recognized.") }
  if(!is.na(gage.type)) { bio.env.dat <- bio.env.dat[grepl(paste0("^", gage.type, "$"), bio.env.dat$type, ignore.case = TRUE), ] }
  
  ##set up env vars ----
  ##- remove NA env var rows
  # na.env <- bio.env.dat %>% filter(if_any(all_of(env.variable.names), is.na))
  bio.env.dat <- bio.env.dat %>% filter(!if_any(all_of(env.variable.names), is.na))

  env_col <- bio.env.dat[, grepl(paste(paste0("^", env.variable.names, "$"), collapse = "|"), colnames(bio.env.dat))]
  if(any(!env.variable.names %in% colnames(env_col))) {
    warning(paste("The following environmental variables were not in the dataframe:", paste(env.variable.names[!env.variable.names %in% colnames(env_col)], collapse = ", ")))
  }
  ## remove constant env vars
  col_rem <- sapply(env_col, function(col) all(col == col[1]))
  env_col <- env_col[, !col_rem]
  
  ##pull out species ----
  spp_col <- bio.env.dat[, grepl(paste(species.list, collapse = "|"), colnames(bio.env.dat))]
  ##- remove species w/ < 10 records
  if(bio.type == "taxonomic") {
    spp_col <- spp_col[, colSums(spp_col) >= 10]
    # ##- remove species w/ < 5% of max collection records?? ##not sure they did this in the end...
    # spp_col <- spp_col[, colSums(spp_col) >= max(colSums(spp_col))*0.05]
    ##- remove species present at every site
    spp_col <- spp_col[, colSums(spp_col) != nrow(spp_col)]
  }
  
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
  
  # gf.out$call$nbin <- nbin
  # gf.out$call$compact <- compact
  
  return(gf.out)
}

#run multiple flow types ----
dir.create(paste0(PATH, "/10_GFOutput/", gsub("-", "_", Sys.Date())))

#break out flows
flow.opt <- c("Int", "RO", "GW")

# for(i in seq_along(site.list)) {
  
  bio <- "fish"
  # if(i == 2) { bio <- "fish"} else { bio <- "bugs" } #set taxa
  
  # comb.dat <- left_join(site.list[[i]], env.dat %>% select(-flw_type, -COMID), by = c("gage_no_15yr" = "site_no"))
  # spp <- names(site.list[[i]] %>% select(-c(lat, long, COMID, contains("flw"), site_id, contains("gage"))))
  
  for(f in flow.opt) {
    
    cat("\nRunning", bio, ":", f, "flow...\n")
   
    file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_traits_", f, "_all.rds")
    if(file.exists(file.name)) { cat("Model exists, next.\n"); next }
    
    gf.out <- NULL
  system.time(
    gf.out <- gf_sub(bio.env.dat = site.list[[3]], #biological data + env data dataframe
                     bio.type = "trait", #if biological vs. trait data, either "trait" or "taxonomic"
                     species.names = species.list, #species of interest (in case less than the full dataframe)
                     env.variable.names = var.names, #variables of interest; previously pulled out HIT, LULC, and full run
                     flow.type = f, #flow type to run, set as NA if all
                     gage.type = NA, #gage type to run, set as NA if all
                     gf.classification = FALSE, #classification or regression? think regression when using trait abundance
                     #gradientForest options:
                     ntree = 999,
                     transform = NULL,
                     compact = TRUE,
                     trace = TRUE,  
                     nbin = 201,
                     corr.threshold = 0.5)
  )
    
    saveRDS(gf.out, file = file.name) #saving bc it would be massive to store all of them in memory

    Sys.sleep(1)
    
  }
# }



#run several groups of models at once: ----
#models to run: 
#fish + bug (2)
# all[[1]], all[[2]]
#all flows, int, ro, gw
# flow.opt <- c("Int", "RO", "GW")
# #all gages, ref, non-ref
# gage.opt <- c(NA, "ref", "non-ref")
# #all vars, HIT, LULC, H_ALT
# var.opt <- list(vars.all$var, 
#                 vars.all[vars.all$var_type == "HIT", ]$var, 
#                 vars.all[vars.all$var_type == "LULC", ]$var, 
#                 vars.all[vars.all$var_type == "H_ALT", ]$var)

#run models
# dir.create(paste0(PATH, "/10_GFOutput/", gsub("-", "_", Sys.Date()))) #date directory for future comparisons
# 
# for(i in 1:length(all)) {
#   for(f in flow.opt) {
#     for(g in gage.opt) {
#       for(v in 1:length(var.opt)) {
#         
#         if(v == 1) { var <- "all" }
#         if(v == 2) { var <- "HIT" }
#         if(v == 3) { var <- "LULC" }
#         if(v == 4) { var <- "HyALT" }
#         if(i == 1) { bio <- "fish"} else { bio <- "bugs" }
#         if(is.na(g)) { ga <- "all" } else { ga <- g }
#         
#         cat("\nRunning", bio, ":", f, "flow,", ga, "gage,", var, "variable type...\n")
#         # file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_", f, "_", ga, "gage_", var, ".rds")
#         file.name <- paste0(PATH,  "/10_GFOutput/", gsub("-", "_", Sys.Date()), "/gf_", bio, "_", f, "_", ga, "gage_", var, ".rds") #forgot x y vars
#         if(file.exists(file.name)) { cat("Model exists, next.\n"); next }
#         
#         gf.out <- NULL
# 
#         gf.out <- gf_sub(bio.env.dat = all[[i]], #biological data + env data dataframe
#                          species.list = spp_list$col_new, #species of interest (in case less than the full dataframe)
#                          env.variable.names = var.opt[[v]], #variables of interest; previously pulled out HIT, LULC, and full run
#                          flow.type = f, #flow type to run, set as NA if all
#                          gage.type = g #gage type to run, set as NA if all
#         )
# 
#         saveRDS(gf.out, file = file.name) #saving bc it would be massive to store all of them in memory
# 
#         Sys.sleep(1)
# 
#       }
#     }
#   }
# }

#extra
# gf.out$call$nbin <- 201
# plot(gf.out)
# most_important <- names(importance(gf.out))[1:25]
# 
# par(mgp = c(2, 0.75, 0))
# plot(gf.out, plot.type = "S", imp.vars = most_important,
#      leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
#      cex.lab = 0.7, line.ylab = 0.9,
#      par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
# 
# plot(gf.out, plot.type = "C", imp.vars = most_important,
#      show.overall = F, legend = T, leg.posn = "topleft",
#      leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4,
#      cex.axis = 0.6, line.ylab = 0.9, 
#      par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
# 
# plot(gf.out, plot.type = "P", show.names = T, horizontal = F,
#      cex.axis = 1, cex.labels = 0.7, line = 2.5)
# 
# plot(gf.out, plot.type = "C", imp.vars = most_important,
#      show.species = F, common.scale = T, cex.axis = 0.6,
#      cex.lab = 0.7, line.ylab = 0.9, 
#      par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
# 
# plot(bio.env.dat$WSAREASQKM, bio.env.dat$longevity)
# plot(bio.env.dat$KFFACTWS, bio.env.dat$temp_pref_cold)

