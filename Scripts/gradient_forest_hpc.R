### GRADIENT FORESTS ###
#code altered to run on HPC, with separate jobs for each bio type x flow type x env type
#10_gradient_forests.R is for local runs looping through options (VERY slow for trait abundance)

suppressPackageStartupMessages({
  library(tidyverse)
  library(gradientForest)
})

options(readr.show_col_types = FALSE)

PATH <- getwd()
message(paste0("\n\nworking directory: ", PATH, "\n\n"))

#read + prep data ----
#for adding watershed info to bio data
huc_nhd_df <- read_csv(paste0(PATH, "/nhd_huc_intersection_info.csv"), col_types = cols(prop_in_huc = col_number(), .default = col_character())) 

#read in bio data (file specified in shell script)
bio.file <- "/BIO_FILE_HERE"
bio.sites <- read_csv(paste0(PATH, bio.file), col_types = cols(COMID = col_character()))

#fix names
names(bio.sites) <- gsub("-", "_", names(bio.sites))

if(any(names(bio.sites) == "gage_no_15yr")) {
  names(bio.sites)[names(bio.sites) == "gage_no_15yr"] <- "site_no"
}

#add in watershed info for later joining
bio.sites <- left_join(bio.sites, huc_nhd_df[, names(huc_nhd_df) %in% c("comid", "huc_id")], by = c("COMID" = "comid")) %>%
  rename(huc12 = huc_id)

species.list <- names(bio.sites)[!names(bio.sites) %in% c("lat", "long", "COMID", "flw_type", "gage_no_15yr", "site_no", "dist2gage_m_15yr", "dist2strm_m_flw", "site_id", "huc12")] 

#env dat
env.info <- read_csv(paste0(PATH, "/environmental_variable_info.csv"))
env.info <- env.info[env.info$VAR_SELECTION_COL == "yes" & !is.na(env.info$VAR_SELECTION_COL), ] #only include previously selected variables

env.type <- gsub("include_", "", "VAR_SELECTION_COL")

var.names <- env.info$variable_col_name

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

env.list <- lapply(var.source[var.source != "taxa dataset"], function(x) {
  x <- gsub("/02_EnvDat", "", x)
  file.name <- paste0(PATH, x)
  env.file <- read_csv(file.name)
  env.file <- env.file %>% mutate(across(matches("COMID|site_no|huc12"), as.character))
  env.file <- env.file[, names(env.file) %in% c(var.names, "huc12", "season", "nlcd_class", "COMID", "site_no")]
  
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
  
  return(env.file)
})

#combine bio + env dat
for(i in seq_along(env.list)) {
  bio.sites <- left_join(bio.sites, env.list[[i]])
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
if(grepl("trait", bio.file, ignore.case = TRUE)) {
  set.bio.type <- "trait"
  
  #set trait file naming option
  if(grepl("cont2cat", bio.file)) { bio.name <- "cont2cat" }
  if(grepl("cont2mn", bio.file)) { bio.name <- "cont2mn" }
  if(grepl("clust", bio.file)) { bio.name <- "clust" }
  
} else {
  set.bio.type <- "taxonomic"
  bio.name <- "pa"
}

#regression if abundance, else classification
if(grepl("abundance", bio.file, ignore.case = TRUE)) {
  set.gf.class <- FALSE
} else {
  set.gf.class <- TRUE
}

#make function for repeated runs
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
  
  gf.out$call$nbin <- nbin
  gf.out$call$compact <- compact
  
  return(gf.out)
}

#run GF ----
# dir.create(paste0(PATH, "/", gsub("-", "_", Sys.Date()))) #set save folder

message("\nRunning ", bio, " x ", set.bio.type, " x ", bio.name, " x ", flow.opt, " flow x ", env.type, "...\n")
   
file.name <- paste0(PATH,  "/output/gf_", bio, "_", set.bio.type, "_", bio.name, "_", flow.opt, "_", env.type, ".rds")

# message("included variables are: ", paste0(var.names, sep = ", "))
# message("included species units are: ", paste0(species.list, sep = ", "))
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
