## Clean + compile occurrence data for fish and inverts

library(tidyverse); library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

dat_path <- paste0(PATH, "/01_BioDat/raw_occurrence_data")

#ADEQ data ----
file.list <- list.files(dat_path, pattern = "^ADEQ.*\\.csv$", full.names = T)
shp.list <- list.files(dat_path, pattern = "ADEQ.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

adeq.l <- lapply(file.list, read_csv)

lapply(adeq.l, problems) #find load errors
# adeq.l[[4]][1190, 4] ##fixed

#remove duplicate columns
adeq.l <- lapply(adeq.l, function(df){
  df <- df[!duplicated(t(df))]
  return(df)
})

## clean bugs ----
bug.l <- adeq.l[grepl("bug", file.list)]
bug.shp <- shp.list[grepl("bug", shp.list)]

unique(unlist(lapply(bug.l, colnames))) #look at all columns to find potential selections

bug.l <- lapply(bug.l, function(df){
  df <- df %>% 
    select(contains(c("stationid", "date", "lat", "long", "utm",
                      "insect", "finalid", "first", "individual"))) %>%
    rename_with(~gsub(".*", "date", .x), contains("startdate")) %>%
    rename_with(~gsub(".*", "taxa_order", .x), contains("orde")) %>%
    rename_with(~gsub(".*", "taxa_name", .x), contains(c("finalid", "first"))) %>%
    rename_with(~gsub(".*", "taxa_count", .x), contains("individ")) %>%
    mutate(date = ifelse(grepl(" ", date), as.Date(date, format = "%m/%d/%Y %H:%M:%S"), mdy(date)),
           date = as.Date(date),
           source = "ADEQ") %>%
    mutate(across(matches("taxa_name"), ~ gsub("\\s.*", "", .x)))
  return(df)
})

bug.l <- bind_rows(bug.l)

bug.shp <- lapply(bug.shp, function(shp){
  df <- read_sf(dat_path, layer = shp)
  df <- df %>%
    st_drop_geometry() %>%
    select(contains(c("stationid", "date", "lat", "long", "utm",
                      "insect", "finalid", "first", "individual")))
  df <- df[!duplicated(t(df))] #remove potential duplicate date columns
  df <- df %>%
    rename_with(~gsub(".*", "date", .x), contains("startdate")) %>%
    rename_with(~gsub(".*", "taxa_order", .x), contains("orde")) %>%
    rename_with(~gsub(".*", "taxa_name", .x), contains(c("finalid", "first"))) %>%
    rename_with(~gsub(".*", "taxa_count", .x), contains("individ")) %>%
    mutate(date = as.Date(date),
           source = "ADEQ") %>%
    mutate(across(matches("taxa_name"), ~ gsub("\\s.*", "", .x)))
  return(df)
})

bug.shp <- bind_rows(bug.shp)

bug.l <- bind_rows(bug.l, bug.shp)

#get station locations to add to NA coordinate rows
##note: "ADEQ3J-21" has multiple taxa sampled, but no lat long included and isn't in the ADEQ_Stations_23.csv file
stid <- bug.l %>% filter(!is.na(Lat)) %>% select(StationID, Lat, Long) %>% distinct()

na.coord <- bug.l %>% filter(is.na(Lat)) %>% select(-Lat, -Long)

na.coord <- left_join(na.coord, stid, by = "StationID")

bug.l <- bind_rows(bug.l %>% filter(!is.na(Lat)), na.coord) #add na coordinate rows back in with coordinates

##do a similar thing with counts, if na for count, see if there's a match (date, species, station) without count
na.ct <- bug.l %>% filter(is.na(taxa_count)) %>% select(-taxa_count)
ctid <- bug.l %>% filter(!is.na(taxa_count)) %>% select(StationID, taxa_name, date, taxa_count) %>% distinct()
na.ct <- na.ct %>%
  left_join(ctid, by = c("StationID", "taxa_name", "date"), multiple = "all") %>%
  distinct() %>%
  group_by(StationID, taxa_name, date) %>%
  mutate(taxa_count = sum(taxa_count)) #for some reason some surveys have multiple counts, going to combine for now

bug.l <- bind_rows(bug.l %>% filter(!is.na(taxa_count)), na.ct)

bug.l <- distinct(bug.l) #remove duplicates
bug.l <- filter(bug.l, !is.na(taxa_name)) #remove records with no taxa (file wasn't for that apparently)

#save combined and cleaned data
write_csv(bug.l, paste0(PATH, "/01_BioDat/source_dat_cleaned/bugs_ADEQ_clean_", gsub("-", "", Sys.Date()), ".csv"))

## clean fish ----
fish.l <- adeq.l[grepl("fish", file.list, ignore.case = T)]
fish.shp <- shp.list[grepl("fish", shp.list, ignore.case = T)]

unique(unlist(lapply(fish.l, colnames))) #look at all columns to find potential selections

fish.l <- lapply(fish.l, function(df){
  df <- df %>% 
    select(contains(c("stationid", "date", "lat", "long",
                      "family", "species", "individ")),
           -contains("KeySpecies")) %>%
    mutate(across(contains("date"), ~as.Date(.x, tryFormats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y"))))
  df <- df[!duplicated(t(df))] #remove potential duplicate columns
  
  #handle the date columns
  date_cols <- df %>% select(contains("date"))
  df$date_min <- do.call(pmin, c(date_cols, na.rm = TRUE))
  df$date_max <- do.call(pmax, c(date_cols, na.rm = TRUE))
  df$date_equal <- df$date_min == df$date_max
  df <- df %>% select(-contains(c("CollDate", "TTL")))
  
  #reformat column names etc.
  df <- df %>%
    rename_with(~gsub(".*", "taxa_family", .x), contains("family")) %>%
    rename_with(~gsub(".*", "taxa_name", .x), contains(c("species"))) %>%
    rename_with(~gsub(".*", "taxa_count", .x), contains("individ")) %>%
    mutate(source = "ADEQ") %>%
    mutate(taxa_genus = gsub("\\s.*", "", taxa_name))
  return(df)
})

fish.l <- bind_rows(fish.l)

fish.shp <- lapply(fish.shp, function(shp){
  df <- read_sf(dat_path, layer = shp)
  df <- df %>%
    st_drop_geometry() %>%
    select(contains(c("stationid", "date", "lat", "long",
                      "family", "species", "individ")),
           -contains("KeySpecies")) %>%
    mutate(across(contains("date"), ~as.Date(.x, tryFormats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y"))))
  df <- df[!duplicated(t(df))] #remove potential duplicate columns
  
  #handle the date columns
  date_cols <- df %>% select(contains("date"))
  df$date_min <- do.call(pmin, c(date_cols, na.rm = TRUE))
  df$date_max <- do.call(pmax, c(date_cols, na.rm = TRUE))
  df$date_equal <- df$date_min == df$date_max
  df <- df %>% select(-contains(c("CollDate", "TTL")))
  
  #reformat column names etc.
  df <- df %>%
    rename_with(~gsub(".*", "taxa_family", .x), contains("family")) %>%
    rename_with(~gsub(".*", "taxa_name", .x), contains(c("species"))) %>%
    rename_with(~gsub(".*", "taxa_count", .x), contains("individ")) %>%
    mutate(source = "ADEQ") %>%
    mutate(taxa_genus = gsub("\\s.*", "", taxa_name))
  return(df)
})

fish.shp <- bind_rows(fish.shp)

fish.l <- bind_rows(fish.l, fish.shp)

#get station locations to add to NA coordinate rows
##note: "ADEQ3J-21" has multiple taxa sampled, but no lat long included and isn't in the ADEQ_Stations_23.csv file
stid <- fish.l %>% filter(!is.na(Lat) & Lat != 0 & !is.na(Long) & Long != 0) %>% select(StationID, Lat, Long) %>% distinct()

na.coord <- fish.l %>% filter(is.na(Lat) | Lat == 0 | is.na(Long) | Long == 0) %>% select(-Lat, -Long)

na.coord <- left_join(na.coord, stid, by = "StationID")

fish.l <- bind_rows(fish.l %>% filter(!(is.na(Lat) | Lat == 0 | is.na(Long) | Long == 0)), na.coord) #add na coordinate rows back in with coordinates

##do a similar thing with counts, if na for count, see if there's a match (date, species, station) without count
na.ct <- fish.l %>% filter(is.na(taxa_count)) %>% select(-taxa_count)

ctid <- fish.l %>% filter(!is.na(taxa_count)) %>% select(StationID, taxa_name, contains("date"), taxa_count) %>% distinct()

na.ct <- na.ct %>%
  left_join(ctid, by = c("StationID", "taxa_name", "date_min", "date_max", "date_equal"), multiple = "all") %>%
  distinct() %>%
  group_by(StationID, taxa_name, date_min, date_max, date_equal) %>%
  mutate(taxa_count = sum(taxa_count)) #for some reason some surveys have multiple counts, going to combine for now

fish.l <- bind_rows(fish.l %>% filter(!is.na(taxa_count)), na.ct)

fish.l <- distinct(fish.l) #remove duplicates
fish.l <- filter(fish.l, !is.na(taxa_name)) #remove records with no taxa (file wasn't for that apparently)

#save combined and cleaned data
write_csv(fish.l, paste0(PATH, "/01_BioDat/source_dat_cleaned/fish_ADEQ_clean_", gsub("-", "", Sys.Date()), ".csv"))


#NAQWA data ----
rm(list = ls()[!ls() %in% c("PATH", "dat_path")])

file.list <- list.files(dat_path, pattern = "USGS.*\\.csv$", full.names = T)
shp.list <- list.files(dat_path, pattern = "USGS.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

usgs.l <- lapply(file.list, read_csv)

lapply(usgs.l, problems) #find load errors, looking at it, none of the problematic columns will be kept

#remove duplicate columns
usgs.l <- lapply(usgs.l, function(df){
  df <- df[!duplicated(t(df))]
  return(df)
})

## clean bugs ----
bug.l <- usgs.l[grepl("bug|insect|invert", file.list, ignore.case = T)]
bug.shp <- shp.list[grepl("bug|insect|invert", shp.list, ignore.case = T)] #no bug usgs shapefiles

unique(unlist(lapply(bug.l, colnames))) #look at all columns to find potential selections

bug.l <- lapply(bug.l, function(df){
  # df <- bug.l[[1]] #testing...
  df <- df %>% 
    select(contains(c("sitenumber", "long", "lat", "collectiondate", "abundance", "coordinate")), 
           matches(c("^PublishedTaxonName$", "^Family$", "^Order$", "^Genus$"))) %>%
    mutate(SiteNumber = as.character(SiteNumber),
           SiteNumber = ifelse(nchar(SiteNumber) < 8, paste0("0", SiteNumber), SiteNumber))
  
if(any(grepl("long", colnames(df), ignore.case = TRUE))) {
    #fix coordinates, get rid of double columns
    ll <- df %>% 
      select(SiteNumber, contains(c("long", "lat", "coordinate"))) %>% 
      distinct()
    ll_cols <- c(
        unique(colnames(ll %>% select(contains("lat")))[max.col(is.na(ll %>% select(contains("lat"))), "first")]),
        unique(colnames(ll %>% select(contains("long")))[max.col(is.na(ll %>% select(contains("long"))), "first")])
        )
    ll <- ll %>%
      select(-all_of(ll_cols)) %>%
      distinct() %>%
      rename_with(~gsub(".*", "Long", .x), contains("long")) %>%
      rename_with(~gsub(".*", "Lat", .x), contains("lat"))
    
    df <- df %>% 
      select(-contains(c("long", "lat", "coordinate"))) %>%
      left_join(., ll, by = "SiteNumber")
  }
  
  #fix date and col names
  df <- df %>% 
    rename_with(~gsub(".*", "date", .x), contains("date")) %>%
    rename_with(~gsub(".*", "taxa_family", .x), contains("family")) %>%
    rename_with(~gsub(".*", "taxa_order", .x), contains("orde")) %>%
    rename_with(~gsub(".*", "taxa_genus", .x), contains("genus")) %>%
    rename_with(~gsub(".*", "taxa_name", .x), contains(c("PublishedTaxonName"))) %>%
    rename_with(~gsub(".*", "taxa_count", .x), contains("abund")) %>%
    mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y %H:%M:%S", "%m/%d/%Y")),
           source = "USGS")
  
    return(df)
})

bug.l <- bind_rows(bug.l)

bug.l <- bug.l %>%
  # filter(!is.na(Long) & !is.na(Lat)) %>%  #remove rows without coordinates
  mutate(taxa_count = ifelse(is.na(taxa_count), 0, taxa_count)) %>% #based on abundance range, NAs likely 0
  filter(!if_all(contains(c("taxa_name", "taxa_family", "taxa_order", "taxa_genus")), ~is.na(.))) %>%
  distinct()

#fix NA coords
stid <- bug.l %>% filter(!is.na(Lat)) %>% select(SiteNumber, Lat, Long, CoordinateDatum) %>% distinct()

na.coord <- bug.l %>% filter(is.na(Lat)) %>% select(-Lat, -Long, -CoordinateDatum)

na.coord <- left_join(na.coord, stid, by = "SiteNumber")

bug.l <- bind_rows(bug.l %>% filter(!is.na(Lat)), na.coord) #add na coordinate rows back in with coordinates


#save combined and cleaned data
write_csv(bug.l, paste0(PATH, "/01_BioDat/source_dat_cleaned/bugs_USGS_clean_", gsub("-", "", Sys.Date()), ".csv"))



# differences in gage data ----

Gages_and_years_with_adjustments <- read_excel("~/Downloads/Gages and years with adjustments.xlsx", 
                                               col_types = c("numeric", "numeric", "text", 
                                                             "numeric", "text", "text", "numeric", 
                                                             "text", "numeric", "numeric", "numeric", 
                                                             "numeric", rep("numeric", 181)), na = c("NA", "NaN", "Inf"))
gage_year_HDI_202308XX <- read_excel("02_EnvDat/gage_year_HDI_202308XX.xlsx",
                                     col_types = c("numeric", "numeric", "text", 
                                                   "numeric", "text", "text", "numeric", 
                                                   "text", "numeric", "numeric", "numeric", 
                                                   "numeric", rep("numeric", 181)), na = c("NA", "NaN", "Inf"))

tmp <- bind_rows(
  setdiff(Gages_and_years_with_adjustments, gage_year_HDI_202308XX),
  setdiff(gage_year_HDI_202308XX, Gages_and_years_with_adjustments)
)
                                                             