## Clean + compile occurrence data for fish and inverts

library(tidyverse); library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

dat_path <- paste0(PATH, "/01_BioDat/raw_occurrence_data")

source(paste0(PATH, "/Scripts/XX_occ_clean_func.R")) ##requires tidyverse//dplyr

#ADEQ data ----
## convert shp to csv
# shp.list <- list.files(dat_path, pattern = "ADEQ.*\\.shp", full.names = F) %>%
#   sub("\\.shp.*$", "", .) %>%
#   unique()
# 
# lapply(shp.list, function(x) {
#   st_read(dsn = dat_path, layer = x) %>%
#     st_drop_geometry() %>%
#     write_csv(., file = paste0(dat_path, "/", x, ".csv"))
# })

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "^ADEQ.*\\.csv$", full.names = T)

# adeq <- lapply(file.list, function(x) { fix_clean_occ(x, set.source = "ADEQ", save.it = NA) })
adeq <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  adeq[[i]] <- fix_clean_occ(file.list[[i]], species.id.col = c("First", "family", "orde", "genus", "species", "FinalID", "^PublishedTaxonName$"),
                             set.source = "ADEQ", save.it = NA)
  
  flush.console()
}

adeq <- bind_rows(adeq)

#fix missing coordinates or incorrect coordinates
adeq <- adeq %>% 
  mutate(lat = ifelse(lat == 0 | long == 0, NA, lat), #change 0 coordinates to NA to try and fix next
         long = ifelse(lat == 0 | long == 0, NA, long), 
         lax = ifelse(is.na(lat), NA, lax), #correct difference from expected columns
         lox = ifelse(is.na(long), NA, lox),
         across(c(lax, lox), ~ifelse(.x == 1, NA, .x))) %>% #if 1, no change needs to occur for coordinates
         distinct() #get rid of duplicates across all columns

unique(adeq$lax) #no lat coordinates to fix
unique(adeq$lox) #try to fix with replacements first

adeq <- adeq %>% #fix extremely diff coordinates (wrong decimal place)
  mutate(long = case_when(lox == 10 ~ (long / 10),
                          lox == 911255 ~ (long / 1000000),
                          T ~ long),
         lox = ifelse(long > 900, "fix", NA)) #update for tracking

stid <- adeq %>% filter(!is.na(lat) & !is.na(long)) %>% select(site_id, lat, long) %>% distinct() #get unique sites w/ coords

na.coord <- adeq %>% filter(is.na(lat) & is.na(long)) %>% select(-lat, -long)

na.coord <- left_join(na.coord, stid, by = "site_id") #add coordinates back in

adeq2 <- bind_rows(adeq %>% filter(!is.na(lat) & !is.na(long)), na.coord) #add na coordinate rows back in with coordinates

nrow(adeq2) == nrow(adeq) #make sure no rows were lost

adeq <- adeq2 %>% 
  select(-lox, -lax) %>% #remove since corrected
  mutate(date = case_when(date_equal == TRUE ~ date_max, #add date for single date surveys
                          date_equal == FALSE ~ date,
                          is.na(date_equal) ~ date),
         across(c(date_max, date_min), ~ case_when(date_equal == TRUE ~ .x == NA,
                                                   T ~ .x)),
         date_equal = ifelse(date_equal == TRUE, NA, date_equal)) %>%
  filter(!complete.cases(order, taxa_name, family, species)) %>% #remove rows with no species id'ed
  filter(!is.na(lat) & !is.na(long)) %>%
  filter(!if_all(c("order", "taxa_name", "family", "species"), is.na)) %>%
  mutate(taxa_name = case_when(bio_type == "bug" ~ sub("\\s.*", "", taxa_name),
                               bio_type == "fish" ~ species)) %>%
  distinct() #check for duplicates again just in case

#save combined and cleaned data
write_csv(adeq, paste0(PATH, "/01_BioDat/source_dat_cleaned/ADEQ_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
# df <- read_csv(paste0(PATH, "/01_BioDat/source_dat_cleaned/ADEQ_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(adeq, adeq2, na.coord, stid)

#NAQWA data ----
#convert shp to csv
# shp.list <- list.files(dat_path, pattern = "USGS.*\\.shp|20171018.*\\.shp", full.names = F) %>%
#   sub("\\.shp.*$", "", .) %>%
#   unique()
# 
# lapply(shp.list, function(x) {
#   st_read(dsn = dat_path, layer = x) %>%
#     st_drop_geometry() %>%
#     write_csv(., file = paste0(dat_path, "/", x, ".csv"))
# })

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "USGS.*\\.csv$|20171018.*\\.csv$", full.names = T)

usgs <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", i, ".", file.list[[i]], "\n")
  flush.console()
  
  usgs[[i]] <- fix_clean_occ(file.list[[i]], ignore.parse.err = TRUE, set.source = "USGS", save.it = NA)
}

usgs <- bind_rows(usgs)

#fix missing coordinates or incorrect coordinates
# library(dataRetrieval)
unique(nchar(usgs$site_id))
usgs <- usgs %>%
  mutate(site_id = case_when(nchar(site_id) == 15 ~ site_id,
                             !startsWith(site_id, "0") ~ paste0("0", site_id),
                             T ~ site_id)) %>%
  distinct() #get rid of duplicates across all columns

stid <- usgs %>% select(site_id) %>% distinct() #get unique sites

stid.info <- dataRetrieval::whatNWISsites(sites = stid$site_id)
stid.info <- stid.info %>%
  rename(site_id = site_no,
         lat = dec_lat_va,
         long = dec_long_va) %>%
  select(site_id, lat, long)

#figure out missing sites
missing.site <- stid %>% filter(!site_id %in% stid.info$site_no)
usgs %>% filter(site_id %in% missing.site$site_id) %>% select(site_id, long, lat) %>% distinct() %>% View()
stid.info <- rbind(stid.info, usgs %>% 
  filter(site_id %in% missing.site$site_id) %>% 
  select(site_id, long, lat) %>% 
  distinct() %>%
  group_by(site_id) %>%
  summarise(lat = mean(lat, na.rm = T),
            long = mean(long, na.rm = T)))

usgs <- usgs %>%
  select(-lat, -long) %>%
  left_join(., stid.info) %>% #add back in coordinates
  distinct() #remove duplicates now that lat/long is corrected

usgs <- usgs %>% 
  filter(!is.na(lat) & !is.na(long)) %>%
  filter(!if_all(c("taxa_name", "order", "family", "genus", "species"), is.na)) %>%
  filter(HybridFlag == "N" | is.na(HybridFlag)) %>%
  group_by(across(-bio_type)) %>%
  slice_head() #check for duplicates again just in case, remove duplicates because of bio_type diff

#save combined and cleaned data
write_csv(usgs, paste0(PATH, "/01_BioDat/source_dat_cleaned/usgs_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(usgs, missing.site, stid, stid.info, shp.list)

#CARLISLE ----
#clean data + bind + save
file.list <- list.files(dat_path, pattern = "carlisle.*\\.csv$", full.names = T)

carl <- list()
carl[[1]] <- fix_clean_occ(file.list[[1]], 
                           site.id.col = "(?<!\\.)\\bID\\b",
                           species.id.col = "taxon",
                           count.id.col = "total",
                           set.source = "Carlisle", save.it = NA)
carl[[2]] <- fix_clean_occ(file.list[[2]], 
                           site.id.col = "ID",
                           species.id.col = "taxon",
                           count.id.col = "total",
                           set.source = "Carlisle", save.it = NA)
# b <- fix_clean_occ(file.list[[3]], ignore.parse.err = FALSE,
#                            site.id.col = "UID",
#                            species.id.col = "taxon",
#                            count.id.col = "total",
#                            set.source = "Carlisle", save.it = NA)
carl[[3]] <- new_df
carl[[4]] <- fix_clean_occ(file.list[[4]], 
                           site.id.col = "ID",
                           species.id.col = "taxon",
                           count.id.col = "total",
                           set.source = "Carlisle", save.it = NA)

carl <- bind_rows(carl) %>% 
  mutate(taxa_name = str_to_sentence(taxa_name))

#fix missing coordinates or incorrect coordinates
unique(carl$lox); unique(carl$lax) #barely out, ignore
carl <- carl %>% select(-lox, -lax)

stid <- carl %>% select(site_id, lat, long, Year) %>% filter(!is.na(lat) & !is.na(long)) %>% distinct() #get unique sites
stid[duplicated(stid$site_id),]

na.coord <- carl %>% filter(is.na(lat) & is.na(long)) %>% select(-lat, -long, -Year)
na.coord <- left_join(na.coord, stid, multiple = "all") %>%
  group_by(across(-Year)) %>%
  summarise(date_min = min(Year, na.rm = TRUE),
            date_max = max(Year, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(contains("date"), ~ifelse(. == Inf | . == -Inf, NA, .)),
         date_equal = date_min == date_max) %>%
  mutate(date = case_when(date_equal == TRUE ~ date_max, #add date for single date surveys
                          date_equal == FALSE ~ NA,
                          is.na(date_equal) ~ NA),
         across(c(date_max, date_min), ~ case_when(date_equal == TRUE ~ .x == NA,
                                                   T ~ .x)),
         date_equal = ifelse(date_equal == TRUE, NA, date_equal))

st.na <- na.coord %>% filter(is.na(lat) & !is.na(taxa_name)) %>% select(site_id) %>% distinct() 
#I can't find any files with coordinates for these sites, so I'm just going to remove them

carl2 <- rbind(na.coord, carl %>%
                 filter(!is.na(lat) & !is.na(long)) %>%
                 mutate(date_min = NA, date_max = NA, date_equal = NA, date = Year) %>% select(-Year))
if(nrow(carl2) == nrow(carl)) {carl <- carl2} else {cat("not equal.\n")}

carl <- carl %>% 
  filter(!is.na(lat) & !is.na(long)) %>%
  filter(!is.na(taxa_name)) %>%
  group_by(across(-bio_type)) %>%
  slice_head() #check for duplicates again just in case, remove duplicates because of bio_type diff

#save combined and cleaned data
write_csv(carl, paste0(PATH, "/01_BioDat/source_dat_cleaned/carlisle_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))


#MATTHEWS ----
file.list <- list.files(dat_path, pattern = "Matthews.*\\.csv$", full.names = T)
# mat <- fix_clean_occ(file.list, site.id.col = "Stream_Name")
mat <- read_csv("01_BioDat/raw_occurrence_data/Matthews_fish_data.csv") %>%
  rename(site_id = Stream_name, Year = Survey_dates, lat = Lat, long = Lon, species = Species) %>%
  select(-State) %>% 
  mutate(bio_type = "fish") %>%
  distinct()
mat1 <- mat %>% 
  filter(grepl("_", Year)) %>%
  mutate(date_min = sub("_.*$", "", Year),
         date_max = sub("^.*?_", "", Year)) %>% select(-Year) %>%
  mutate(across(contains("date"), ~gsub("\\D", "", .x)),
         across(contains("date"), ~ifelse(.x > 60, paste0("19", .x), paste0("20", .x))),
         date_equal = FALSE) 
mat2 <- mat %>% 
  filter(!grepl("_", Year)) %>%
  rename(date = Year)

mat <- bind_rows(mat1, mat2) %>%
  distinct()

write_csv(mat, paste0(PATH, "/01_BioDat/source_dat_cleaned/matthews_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))

rm(list = ls()[grepl("mat", ls())])

#MO UNKNOWN SOURCE ----
mo <- readxl::read_excel("01_BioDat/raw_occurrence_data/MO_Fish_Community_Data_ALL.xlsx")
write_csv(mo, paste0(dat_path, "/MO_Fish_Community_Data_ALL.csv"))

file.list <- list.files(dat_path, pattern = "^MO.*\\.csv$", full.names = T)
file.list <- file.list[!grepl("USGS|2017", file.list)]

mo <- lapply(file.list[1:2], function(x){
  out <- read_csv(x) %>% 
    rename(site_id = RCH_CODE,
           taxa_name = Taxa,
           genus = Genera)
  out1 <- out %>% 
    #this has coordinates that are WAY off like in south america, so fixing those
    filter(UTM_X < 4493502 & UTM_X > 4025253 & UTM_Y < 843651 & UTM_Y > 314088) %>%
    rename(UTM_X = UTM_Y, UTM_Y = UTM_X)
  out2 <- out %>%
    filter(UTM_X > 285100 & UTM_X < 835900 & UTM_Y > 3988500 & UTM_Y < 4491200)
  out <- bind_rows(out1, out2)
  out <- st_as_sf(out, coords = c("UTM_X", "UTM_Y"), crs = 26915, remove = FALSE) %>%
    st_transform(., crs = 4326)
  
  out <- cbind(out, st_coordinates(out)) %>%
    st_drop_geometry() %>%
    rename(lat = Y, long = X) %>%
    mutate(site_id = as.character(site_id), bio_type = "bug", source = "MO_UNK") %>%
    select(-contains("UTM"), -STATE) %>%
    distinct()
})
mo <- bind_rows(mo) %>% distinct()
  
mo2 <- read_csv(file.list[[3]]) %>%
  select(contains("UTM"), SEG_ID, Date, Scientific, Number) %>%
  rename(site_id = SEG_ID, date = Date, taxa_name = Scientific, taxa_count = Number) %>%
  mutate(taxa_name = str_to_sentence(taxa_name)) %>%
  st_as_sf(., coords = c("UTM_E", "UTM_N"), crs = 26915, remove = FALSE) %>%
  st_transform(., crs = 4326)
mo2 <- cbind(mo2, st_coordinates(mo2)) %>%
  st_drop_geometry() %>%
  rename(lat = Y, long = X) %>%
  mutate(site_id = as.character(site_id), bio_type = "fish", source = "MO_UNK", date = as.Date(date)) %>%
  select(-contains("UTM")) %>%
  distinct()

mo3 <- read.csv(file.list[[4]]) %>%
  select(long, lat, RCH_CODE, Taxa, Date) %>%
  rename(site_id = RCH_CODE, taxa_name = Taxa, date = Date) %>%
  mutate(taxa_name = str_to_sentence(taxa_name),
         site_id = as.character(site_id), bio_type = "bug", source = "MO_UNK",
         date = as.Date(format(strptime(date, format = "%d-%b-%y"), "%Y-%m-%d"))) %>%
  distinct()

mo_un <- bind_rows(mo, mo2, mo3)
mo_un <- mo_un %>%
  filter(!is.na(lat) & !is.na(long)) %>%
  filter(!is.na(taxa_name)) %>%
  distinct(taxa_name, genus, long, lat, date, .keep_all = TRUE) %>%
  mutate(HybridFlag = ifelse(grepl(" x ", taxa_name), "Y", "N")) %>%
  filter(HybridFlag == "N") %>% #remove hybrids
  select(-HybridFlag)

write_csv(mo_un, paste0(PATH, "/01_BioDat/source_dat_cleaned/MO_unknwn_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(list = ls()[grepl("mo|out", ls())])


#MONHP ----
monhp <- st_read(dsn = dat_path, layer = "MONHP") %>%
  select(SNAME, DATEFRST, DATELAST, contains("UTM")) %>%
  st_drop_geometry()
monhp <- st_as_sf(monhp, coords = c("UTME", "UTMN"), crs = 26915) %>%
  st_transform(., crs = 4326)
monhp <- cbind(monhp, st_coordinates(monhp)) %>%
  st_drop_geometry() %>%
  rename(taxa_name = SNAME,
         lat = Y,
         long = X,
         date_min = DATEFRST,
         date_max = DATELAST) %>%
  mutate(date_equal = date_min == date_max) %>%
  distinct() %>%
  mutate(date = case_when(date_equal == TRUE ~ date_max, #add date for single date surveys
                          date_equal == FALSE ~ NA,
                          is.na(date_equal) ~ NA),
         across(c(date_max, date_min), ~ case_when(date_equal == TRUE ~ .x == NA,
                                                   T ~ .x)),
         date_equal = ifelse(date_equal == TRUE, NA, date_equal),
         source = "MONHP",
         bio_type = "unknown")

write_csv(monhp, paste0(PATH, "/01_BioDat/source_dat_cleaned/monhp_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(monhp)

# OK_DEQ ----
## convert shp to csv
shp.list <- list.files(dat_path, pattern = "OK_DEQ.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

lapply(shp.list, function(x) {
  st_read(dsn = dat_path, layer = x) %>%
    st_drop_geometry() %>%
    write_csv(., file = paste0(dat_path, "/", x, ".csv"))
})

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "OK_DEQ.*\\.csv$", full.names = T) %>% unique()

okdeq <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  okdeq[[i]] <- fix_clean_occ(file.list[[i]], 
                             species.id.col = c("genus", "...12"),
                             date.id.col = "Field7",
                             site.id.col = c("sites"),
                             set.source = "OK_DEQ", save.it = NA)
  
  flush.console()
}

okdeq <- bind_rows(okdeq)
okdeq <- okdeq %>%
  select(-lax, -lox) %>% #amount out not a concern
  distinct()

write_csv(okdeq, paste0(PATH, "/01_BioDat/source_dat_cleaned/okdeq_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(okdeq)

#ORWB ----
orwb <- fix_clean_occ(paste0(dat_path, "/OK.bugs.filt_ORWB_23.csv"),
                      species.id.col = c("First"),
                      set.source = "ORWB", save.it = NA)
write_csv(orwb, paste0(PATH, "/01_BioDat/source_dat_cleaned/orwb_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(orwb)


#OKCONS ----
## convert shp to csv
shp.list <- list.files(dat_path, pattern = "OKCons.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

lapply(shp.list, function(x) {
  st_read(dsn = dat_path, layer = x) %>%
    st_drop_geometry() %>%
    write_csv(., file = paste0(dat_path, "/", x, ".csv"))
})

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "OKCons.*\\.csv$", full.names = T) %>% unique()

okcons <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  okcons[[i]] <- fix_clean_occ(file.list[[i]],
                               count.id.col = "Total",
                               site.id.col = "WBID",
                               set.source = "OKCons", save.it = NA)
  
  flush.console()
}

okcons <- bind_rows(okcons) %>% distinct()
# okcons %>% filter(is.na(lat) | is.na(long))

write_csv(okcons, paste0(PATH, "/01_BioDat/source_dat_cleaned/okcons_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(okcons)


#OWSB ----

#DO NEXT!!!!!

#ONHI ----

#QUINN ----

#RAM ----

#REQUEST UNKNOWN SRCE ----

#SMB ----

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
                                                             