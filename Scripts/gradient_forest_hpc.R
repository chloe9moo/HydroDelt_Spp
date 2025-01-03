### GRADIENT FORESTS ###
#code altered to run on HPC, with separate jobs for each bio type x flow type x env type
#10_gradient_forests.R is for local runs looping through options (VERY slow for trait abundance)

suppressPackageStartupMessages({
  library(tidyverse)
  library(gradientForest) ##local
  # library(gradientForest, lib.loc = "/scrfs/storage/voorhees/home/R/x86_64-pc-linux-gnu/4.2")
})

options(readr.show_col_types = FALSE)

#sessionInfo()

PATH <- getwd()
# PATH <- paste0(getwd(), "/working_dir") #for testing, manual runs
message(paste0("\n\nworking directory: ", PATH, "\n\n"))

#read + prep data ----
#for adding watershed info to bio data
huc_nhd_df <- read_csv(paste0(PATH, "/nhd_huc_intersection_info.csv"), col_types = cols(prop_in_huc = col_number(), .default = col_character())) 

#read in bio data (file specified in shell script)
bio.file <- "/BIO_FILE_HERE"
bio.sites <- read_csv(paste0(PATH, bio.file), col_types = cols(COMID = col_character()))

##for testing
#print(head(bio.sites))
#stop("TESTING STOP, FILES LOADED")

#fix names
names(bio.sites) <- gsub("-", "_", names(bio.sites))

if(any(names(bio.sites) == "flw_gage_no")) {
  names(bio.sites)[names(bio.sites) == "flw_gage_no"] <- "site_no"
}

#add in watershed info for later joining
bio.sites <- left_join(bio.sites, huc_nhd_df[, names(huc_nhd_df) %in% c("comid", "huc_id")], by = c("COMID" = "comid")) %>%
  rename(huc12 = huc_id)

species.list <- names(bio.sites)[!names(bio.sites) %in% c("lat", "long", "COMID", "flw_type", "flw_gage_no", "site_no", "dist2gage_m_15yr", "dist2strm_m_flw", "site_id", "huc12")] 

#env dat
env.info <- read_csv(paste0(PATH, "/environmental_variable_info.csv"))
env.info <- env.info[env.info$VAR_SELECTION_COL == "yes" & !is.na(env.info$VAR_SELECTION_COL), ] #only include previously selected variables

env.type <- gsub("include_", "", "VAR_SELECTION_COL")

var.names <- env.info$variable_col_name

if(any(grepl("lat|long", var.names))) { lat_long_var <- TRUE } else { lat_long_var <- FALSE } #for adding lat/long as predictor variable later

#select + adjust variable source needed (where necessary)
var.source <- unique(env.info$variable_source)

if(any(grepl("fish|bug", var.source))) {
  
  tmp <- var.source[grepl("fish|bug", var.source)]
  
  if(any(grepl("fish", bio.file))) {
    tmp <- gsub("\\[fish\\|bug\\]", "fish", tmp)
  } else {
    tmp <-  gsub("\\[fish\\|bug\\]", "bug", tmp)
  }
  
  var.source <- c(var.source[!grepl("fish|bug", var.source)], tmp)
  
  rm(tmp)
}

if(any(grepl("all\\|ssn", var.source))) {
  
  tmp <- var.source[grepl("all\\|ssn", var.source)]
  
  tmp <- gsub("\\[all\\|ssn\\]", "ssn", tmp)
  
  var.source <- c(var.source[!grepl("all\\|ssn", var.source)], tmp)
  
  rm(tmp)
}

#load env variable csvs
env.list <- lapply(var.source[var.source != "taxa dataset"], function(x) {
  x <- gsub("/02_EnvDat", "", x)
  file.name <- paste0(PATH, x)
  env.file <- read_csv(file.name)
  env.file <- env.file %>% mutate(across(matches("COMID|site_no|huc12"), as.character))
  env.file <- env.file[, names(env.file) %in% c(var.names, "huc12", "season", "nlcd_class", "COMID", "site_no", "site_id")]
  
  if(any(names(env.file) %in% c("season", "nlcd_class"))) {
    env.file <- env.file %>%
      pivot_longer(cols = -matches("huc12|season|nlcd_class|COMID|site_no"), names_to = "variable", values_to = "value") %>%
      pivot_wider(names_from = c("variable", matches("nlcd_class|season")), values_from = value)
  }
  
  #remove zero sum or all NA columns (won't be useful)
  env.file <- env.file[, !sapply(env.file, function(x) {
    if(is.numeric(x)) {
      all(is.na(x)) | (sum(x, na.rm = TRUE) == 0)
    } else {
      all(is.na(x))
    }
  })]
  
  env.file <- distinct(env.file)
  
  return(env.file)
})

#update variable names for model input (not always necessary but sometimes)
new.names <- c()
for(i in seq_along(env.list)) {
  new.names <- c(new.names, names(env.list[[i]]))
}

var.names <- new.names[!new.names %in% c("COMID", "huc12", "site_no", "site_id")]

if(lat_long_var) { var.names <- c("lat", "long", var.names) }

#combine bio + env dat
##set dataframe stop check (didn't lose anything during join)
ncol_bio <- length(unique(c(names(bio.sites), var.names)))
nrow_bio <- nrow(bio.sites)

huc.sites <- bio.sites[grepl("\\|", bio.sites$huc12), ] #pull out now for later
bio.sites <- bio.sites[!grepl("\\|", bio.sites$huc12), ]

#join regular sites, with single watershed overlap
for(i in seq_along(env.list)) {
  bio.sites <- left_join(bio.sites, env.list[[i]])
}

##deal with nhd streams that overlap more than 1 huc
if(nrow(huc.sites) != 0) {
  #retain huc combo string for joining back later
  bio <- huc.sites[, !names(huc.sites) %in% c("huc12", "COMID", "site_no")]
  ids <- huc.sites[, names(huc.sites) %in% c("site_id", "huc12")]
  
  huc.sites <- huc.sites %>%
    select(site_id, huc12, COMID, site_no) %>%
    separate_longer_delim(huc12, "|")
  
  for(i in seq_along(env.list)) {
    huc.sites <- left_join(huc.sites, env.list[[i]])
  }
  
  #find average for sites across env vars, add combo HUC string back in, add spp data back in
  huc.sites <- huc.sites %>%
    select(-huc12) %>%
    group_by(site_id, COMID, site_no) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA, .x))) %>%
    left_join(., ids) %>%
    left_join(., bio)
  
  bio.sites <- bind_rows(bio.sites, huc.sites)
}

if(nrow(bio.sites) != nrow_bio | ncol(bio.sites) != ncol_bio) { 
  stop("lost sites or columns during env x bio join step...")
}

#run GF models ----
#set flow subset to run
flow.opt <- "FLOW_TYPE_HERE"

#set which taxa for saving
if(grepl("fish", bio.file, ignore.case = TRUE)) {
  bio <- "fish"
}

if(grepl("bug", bio.file, ignore.case = TRUE)) {
  bio <- "bug"
}

#set taxonomic or trait dataset being run
if(grepl("_trait_", bio.file, ignore.case = TRUE)) {
  set.bio.type <- "trait"
  
  #set trait file naming option
  if(grepl("cont2cat", bio.file)) { bio.name <- "cont2cat" }
  if(grepl("cont2mn", bio.file)) { bio.name <- "cont2mn" }
  if(grepl("clust", bio.file)) { bio.name <- "clust" }
  
}
if(grepl("_wide_", bio.file, ignore.case = TRUE)) {
  set.bio.type <- "taxonomic"
  bio.name <- "pa"
}
if(grepl("div", bio.file, ignore.case = TRUE)) {
  set.bio.type <- "div"
  if(grepl("rich", bio.file)) { bio.name <- "rich" }
  if(grepl("fdisp", bio.file)) { bio.name <- "fdisp" }
}

#regression if abundance, else classification
if(grepl("abundance|div", bio.file, ignore.case = TRUE)) {
  set.gf.class <- FALSE
} else {
  set.gf.class <- TRUE
}

#make function for repeated runs (not really necessary here, but keeping anyway for now)
gf_sub <- function(bio.env.dat = bio.sites, #biological data + env data dataframe
                   bio.type = set.bio.type, #if biological vs. trait data, either "trait" or "taxonomic"
                   species.names = species.list, #species of interest (in case less than the full dataframe)
                   env.variable.names = var.names, #variables of interest; previously pulled out HIT, LULC, and full run
                   flow.type = NA, #flow type to run, set as NA if all
                   gage.type = NA, #gage type to run, set as NA if all
                   gf.classification = set.gf.class, #classification or regression? when to use each?
                   #gradientForest options:
                   ntree = 999,
                   transform = NULL,
                   compact = TRUE,
                   trace = TRUE,  
                   nbin = 201,
                   corr.threshold = 0.5) {
  
  ##subset by flow + reference ----
  #flow
  if(!any(grepl(flow.type, c("Int", "RO", "GW"), ignore.case = TRUE) | is.na(flow.type))) { stop("flow type not recognized.") }
  if(!is.na(flow.type)) { bio.env.dat <- bio.env.dat[grepl(flow.type, bio.env.dat$flw_type, ignore.case = TRUE), ] }
  
  #reference
  if(!any(grepl(gage.type, c("ref", "non-ref"), ignore.case = TRUE) | is.na(gage.type))) { stop("gage type not recognized.") }
  if(!is.na(gage.type)) { bio.env.dat <- bio.env.dat[grepl(paste0("^", gage.type, "$"), bio.env.dat$type, ignore.case = TRUE), ] }
  
  ##set up env vars ----
  ##- remove NA env var rows
  na.env <- bio.env.dat %>% filter(if_any(all_of(env.variable.names), is.na)) #%>% select(any_of(env.variable.names))
  print(apply(na.env %>% select(any_of(env.variable.names)), 2, function(x) sum(is.na(x)))) #check where NAs are coming from (is one var more of a culprit?)
  message("Removing ", nrow(na.env), " sites due to NA values in env. vars...\n")
  
  bio.env.dat <- bio.env.dat %>% filter(!if_any(all_of(env.variable.names), is.na))

  env_col <- bio.env.dat[, grepl(paste(paste0("^", env.variable.names, "$"), collapse = "|"), colnames(bio.env.dat))]
  if(any(!env.variable.names %in% colnames(env_col))) {
    warning(paste("The following environmental variables were not in the dataframe:", paste(env.variable.names[!env.variable.names %in% colnames(env_col)], collapse = ", ")))
  }
  ## remove constant env vars
  col_rem <- sapply(env_col, function(col) all(col == col[1]))
  env_col <- env_col[, !col_rem]
  
  ##pull out species/trait columns ----
  spp_col <- bio.env.dat[, grepl(paste(species.list, collapse = "|"), colnames(bio.env.dat))]
  
  ##remove species present at every site or present at none of the sites ----
  if(bio.type == "taxonomic") {
    spp_col <- spp_col[, colSums(spp_col) != nrow(spp_col)]
    spp_col <- spp_col[, colSums(spp_col) != 0]
    # removed_spp <- species.list[!species.list %in% names(spp_col)]
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
  
  gf.out$call$nbin <- nbin
  gf.out$call$compact <- compact
  
  return(gf.out)
}

#run GF ----
# dir.create(paste0(PATH, "/", gsub("-", "_", Sys.Date()))) #set save folder

message("\nRunning ", bio, " x ", set.bio.type, " x ", bio.name, " x ", flow.opt, " flow x ", env.type, "...\n")
   
file.name <- paste0(PATH,  "/output/gf_", bio, "_", set.bio.type, "_", bio.name, "_", flow.opt, "_", env.type, ".rds")

# message("included variables are: ", paste0(var.names, sep = ", "))
# message("included taxa units are: ", paste0(species.list, sep = ", "))
message("\nStart: ", Sys.time(), "\n")

gf.out <- gf_sub(bio.env.dat = bio.sites, #biological data + env data dataframe
                 bio.type = set.bio.type, #if biological vs. trait data, either "trait" or "taxonomic"
                 species.names = species.list, #species of interest (in case less than the full dataframe)
                 env.variable.names = var.names, #variables of interest; previously pulled out HIT, LULC, and full run
                 flow.type = flow.opt, #flow type to run, set as NA if all
                 gage.type = NA, #gage type to run, set as NA if all
                 gf.classification = set.gf.class, #classification or regression? when to use each?
                 #gradientForest options:
                 ntree = 999,
                 transform = NULL,
                 compact = TRUE,
                 trace = TRUE,
                 nbin = 201,
                 corr.threshold = 0.5)

saveRDS(gf.out, file = file.name) #save model output to look at later

message("\nEnd: ", Sys.time(), "\n")

message("\n\nmodel save location: ", file.name, "\n\n")

