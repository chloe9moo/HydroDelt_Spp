## Clean + compile occurrence data for fish and inverts -- function

# library(tidyverse); library(sf)
# options(readr.show_col_types = FALSE)

# PATH <- getwd()

# dat_path <- paste0(PATH, "/01_BioDat/raw_occurrence_data")

# file.list <- list.files(dat_path, pattern = "^ADEQ.*\\.csv$", full.names = T)
# file.list <- list.files(dat_path, pattern = ".csv", full.names = T)

# occ.file.path <- file.list[[3]]

fix_clean_occ <- function(
  occ.file.path = NA,
  ignore.parse.err = FALSE, #after checking problems with load in, ignore?
  site.id.col = c("StationID", "site_no", "SiteNumber", "Site.ID"), #id columns to search for
  date.id.col = c("date", "collection", "year"),
  species.id.col = c("family", "orde", "genus", "species", "FinalID", "^PublishedTaxonName$"), #species columns to search for
  lat.long.check = c(xmax = -85, xmin = -102, ymax = 41, ymin = 30), #expected region for coordinates (for extreme errors)
  count.id.col = c("individ", "abundance"),
  set.source = sub(".csv", "", sub(".*/", "", occ.file.path)),
  save.it = paste0(PATH, "/01_BioDat/source_dat_cleaned/")
  ) {
  
  df <- read_csv(occ.file.path)
  
  if(nrow(problems(df)) != 0 & ignore.parse.err == FALSE) { 
    cat("Problems loading csv..\n")
    return(problems(df))
  }
  
  df <- df[!duplicated(t(df))] #get rid of duplicate columns regardless of col name
  df$track_id <- seq(1, nrow(df))
  
  new_df <- df[, "track_id"]
  tot_row <- nrow(df) #save for later
  
  #add hybrid column for post function filtering
  if(any(grepl("hybrid", colnames(df), ignore.case = T))) {
    new_df <- cbind(new_df, df[, grepl("hybrid", colnames(df), ignore.case = T)])
  }
  
  #get site column ----
  sid <- df[, grepl(paste0(site.id.col, collapse = "|"), colnames(df), ignore.case = T, perl = TRUE) & colnames(df) != "track_id"] 
  sid <- mutate(sid, across(everything(), as.character))
  if(ncol(sid) == 0) { cat("No site id columns found in dataset. \n") }
  colnames(sid) <- "site_id"
  new_df <- cbind(new_df, sid)
  
  #get date column ----
  dc <- df[, grepl(paste0(date.id.col, collapse = "|"), colnames(df), ignore.case = T)] #get date column
  if(ncol(dc) == 0) { cat("No date columns found in dataset. \n") 
    } else {
    dc <- dc[, !grepl(paste0(c("identification", "verification", "curation", "year", "month", "day", "composite"), collapse = "|"), ignore.case = T, colnames(dc))]
    if(ncol(dc) == 0) { 
      dc <- df[, grepl(paste0(date.id.col, collapse = "|"), colnames(df), ignore.case = T)]
      } else {
        dc <- dc %>% mutate(across(everything(), ~sub(" .*", "", .x))) #get rid of time
    
        err <- try(lapply(dc, as.Date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d", "%Y")), silent = T) #check if its already in the correct format
        if(inherits(err, "try-error")) {
          dc <- as.data.frame(lapply(dc, mdy)) #convert to date
        } else {
          dc <- as.data.frame(lapply(dc, as.Date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d", "%Y")))
        }
        
        date_cols <- colnames(dc)
      
        if(ncol(dc) != 1 & any(duplicated(t(dc)) == FALSE)) { #if there is more than one date column AND they are not dups:
          #get range of different date columns
          dc$date_min <- do.call(pmin, c(dc[, date_cols], na.rm = TRUE))
          dc$date_max <- do.call(pmax, c(dc[, date_cols], na.rm = TRUE))
          dc$date_equal <- dc$date_min == dc$date_max
          
        } else { #if there is more than one column that is dupulicated OR only one column
          dc$date <- dc[, 1] #just take the first column and name it date
        }
        
        dc <- dc %>% select(-all_of(date_cols)) #remove original date columns
      }
    
    new_df <- cbind(new_df, dc) #bind to new data frame
  }
    
  #get species id columns ----
  spp <- df[, grepl(paste0(species.id.col, collapse = "|"), colnames(df), ignore.case = TRUE)]
  if(ncol(spp) == 0) { cat("No species columns found in dataset. \n") 
    } else {
    spp <- spp[, !sapply(spp, is.numeric)] #get rid of obviously incorrect ones
    spp <- spp[, !grepl("sub|super|infra|key", colnames(spp), ignore.case = TRUE)]
    spp <- spp %>%
      rename_with(~gsub(".*", "order", .x), contains("orde")) %>%
      rename_with(~gsub(".*", "family", .x), contains("family")) %>%
      rename_with(~gsub(".*", "genus", .x), contains("genus")) %>%
      rename_with(~gsub(".*", "species", .x), contains("species"))
    colnames(spp)[!grepl("order|family|genus|species", colnames(spp))] <- "taxa_name"
    new_df <- cbind(new_df, spp)
  }

  #get lat long ----
  ll <- df[, grepl("lon|lat|datum", colnames(df), ignore.case = TRUE)]
  if(!any(grepl("lon|lat", colnames(ll), ignore.case = TRUE))) { cat("No coordinate columns found in dataset. \n") 
    } else {
    ll <- ll[, !grepl("time|vertical", colnames(ll), ignore.case = TRUE)]
    
    if(sum(grepl("lon|lat", colnames(df), ignore.case = TRUE)) > 2) {
      ll_cols <- c( #remove mostly NA columns
        unique(colnames(ll %>% select(contains("lat")))[max.col(is.na(ll %>% select(contains("lat"))), "first")]),
        unique(colnames(ll %>% select(contains("lon")))[max.col(is.na(ll %>% select(contains("lon"))), "first")])
      )
      ll <- ll[, !grepl(paste0(ll_cols, collapse = "|"), colnames(ll))]
    }
    
    ll <- ll %>%
      rename_with(~gsub(".*", "long", .x), contains("lon")) %>%
      rename_with(~gsub(".*", "lat", .x), contains("lat")) %>%
      rename_with(~gsub(".*", "datum", .x), contains("datum")) %>%
      mutate(long = abs(long)*-1) #make sure it's all negative (western WGS84 expected)
    
    if(any(ll$long > lat.long.check[["xmax"]] | ll$long < lat.long.check[["xmin"]] | 
       ll$lat > lat.long.check[["ymax"]] | ll$lat < lat.long.check[["ymin"]], na.rm = TRUE)) {
      cat("Coordinates not within expected bounding box.\n")
      ll <- ll %>%
        mutate(lox = case_when(long > lat.long.check[["xmax"]] ~ round(long/lat.long.check[["xmax"]]),
                               long < lat.long.check[["xmin"]] ~ round(long/lat.long.check[["xmin"]]),
                               T ~ 1),
               lax = case_when(lat > lat.long.check[["ymax"]] ~ round(long/lat.long.check[["ymax"]]),
                               lat < lat.long.check[["ymin"]] ~ round(long/lat.long.check[["ymin"]]),
                               T ~ 1))
    }
    
    if(any(colnames(ll) == "datum")) {
      crs <- unique(ll$datum) #get datums stored
      ll$ID <- seq(1, nrow(ll))
      
      ll_2 <- data.frame()
      for(i in 1:length(crs)) {
        xy <- ll[ll$datum == crs[i], ]
        xy <- st_as_sf(xy, coords = c("long", "lat"), crs = crs[i]) %>%
          st_transform(xy, crs = 4326) 
        xy <- cbind(xy, st_coordinates(xy)) %>%
          st_drop_geometry() %>%
          rename(long = X, lat = Y) %>% select(-datum)
        ll_2 <- rbind(ll_2, xy)
      }
      
      ll <- ll %>% select(-long, -lat, -datum) %>% #to retain the correct order, join by ID
        left_join(., ll_2) %>%
        select(-ID)
    }
    
    new_df <- cbind(new_df, ll)
  }

  #get abundance / counts ----
  ct <- df[, grepl(paste0(count.id.col, collapse = "|"), colnames(df), ignore.case = TRUE)]
  ct <- ct[, !grepl("TTL", colnames(ct))]
  if(ncol(ct) == 0) { cat("No count columns found in dataset. \n") 
    } else {
      
      if(ncol(ct) > 1) {
        ct$taxa_count <- rowSums(ct, na.rm = T)
        ct <- select(ct, taxa_count)
      } else {
        colnames(ct) <- "taxa_count"
      }
    
    if(min(ct, na.rm = TRUE) != 0) {
      ct[is.na(ct$taxa_count), ] <- 0
    }
    
    new_df <- cbind(new_df, ct)
  }
  
  #set source and type ----
  file.name <- sub(".csv", "", sub(".*/", "", occ.file.path))
  
  new_df <- new_df %>%
    mutate(bio_type = ifelse(grepl("bug|invert", file.name, ignore.case = TRUE), "bug",
                             ifelse(grepl("fish", file.name, ignore.case = TRUE), "fish", 
                                    "unknown")),
           source = set.source)
  
  #distinct & save ----
  if(nrow(new_df) != nrow(df)) { stop("Something went wrong!\nCleaned df and old df don't have the same number of rows.") }
  new_df <- new_df %>%
    select(-track_id) %>%
    distinct()
  
  if(!is.na(save.it)) {
    write_csv(new_df, paste0(save.it, "/", file.name, "_cleaned_", gsub("-", "", Sys.Date()), ".csv"))
  }
  
  return(new_df)
}

# fix_clean_occ(file.list[[1]], set.source = "Carlisle", save.it = NA)
# rm(ct, dc, df, err, ll, new_df, sid, spp, count.id.col, date.id.col, lat.long.check, site.id.col, species.id.col, tot_row)




