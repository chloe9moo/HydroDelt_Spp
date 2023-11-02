## Clean collected occurrence data for fish and inverts

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

#USGS data ----
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

##NOTE! MO_bugs_comb* FILES ARE NOW EXCLUDED. THERE ARE TOO MANY ISSUES IN THE DATASETS
# mo <- lapply(file.list[1:2], function(x){
#   out <- read_csv(x) %>% 
#     rename(site_id = RCH_CODE,
#            taxa_name = Taxa,
#            genus = Genera)
#   out1 <- out %>% 
#     #this has coordinates that are WAY off like in south america, so fixing those
#     filter(UTM_X < 4493502 & UTM_X > 4025253 & UTM_Y < 843651 & UTM_Y > 314088) %>%
#     rename(UTM_X = UTM_Y, UTM_Y = UTM_X)
#   out2 <- out %>%
#     filter(UTM_X > 285100 & UTM_X < 835900 & UTM_Y > 3988500 & UTM_Y < 4491200)
#   out <- bind_rows(out1, out2)
#   out <- st_as_sf(out, coords = c("UTM_X", "UTM_Y"), crs = 26915, remove = FALSE) %>%
#     st_transform(., crs = 4326)
#   
#   out <- cbind(out, st_coordinates(out)) %>%
#     st_drop_geometry() %>%
#     rename(lat = Y, long = X) %>%
#     mutate(site_id = as.character(site_id), bio_type = "bug", source = "MO_UNK") %>%
#     select(-contains("UTM"), -STATE) %>%
#     distinct()
# })
# mo <- bind_rows(mo) %>% distinct()
  
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
  rename(site_id = RCH_CODE, taxa_name = Taxa, date = Date, long = lat, lat = long) %>%
  mutate(taxa_name = str_to_sentence(taxa_name),
         site_id = as.character(site_id), bio_type = "bug", source = "MO_UNK",
         date = as.Date(format(strptime(date, format = "%d-%b-%y"), "%Y-%m-%d"))) %>%
  distinct()

mo_un <- bind_rows(mo2, mo3)
mo_un <- mo_un %>%
  filter(!is.na(lat) & !is.na(long)) %>%
  filter(!is.na(taxa_name)) %>%
  distinct(taxa_name, long, lat, date, .keep_all = TRUE) %>%
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
owsb <- read_csv(paste0(dat_path, "/OKbugs_OWSB.csv")) %>%
  rename(date = 'Activity Start Date',
         lat = 'Monitoring Location Latitude',
         long = 'Monitoring Location Longitude',
         taxa_name = 'Taxonomic Name') %>%
  mutate(date = mdy(date),
         bio_type = "bug",
         source = "OWSB",
         taxa_name = ifelse(grepl("retire", taxa_name, ignore.case = TRUE), sub(".*use ", "", taxa_name), taxa_name)) %>%
  distinct()

write_csv(owsb, paste0(PATH, "/01_BioDat/source_dat_cleaned/owsb_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(owsb)


#ONHI ----
## convert shp to csv
shp.list <- list.files(dat_path, pattern = "ONHI.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

lapply(shp.list, function(x) {
  st_read(dsn = dat_path, layer = x) %>%
    st_drop_geometry() %>%
    write_csv(., file = paste0(dat_path, "/", x, ".csv"))
})

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "ONHI.*\\.csv$", full.names = T) %>% unique()

onhi <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  onhi[[i]] <- fix_clean_occ(file.list[[i]],
                             species.id.col = "sname",
                             count.id.col = "occurrence", #commented out the min to 0 part of the function to make this work (since not numeric)
                             site.id.col = c("location", "locality"),
                             set.source = "ONHI", save.it = NA)
  
  flush.console()
}

onhi <- bind_rows(onhi) %>% distinct()

unique(onhi$lox); unique(onhi$lax) #amount out of box not a concern
onhi <- onhi %>% select(-lox, -lax)

#fix count
onhi <- onhi %>%
  mutate(date = as.character(date),
         date = case_when(grepl("date to year only", taxa_count, ignore.case = TRUE) ~ as.character(sub(".*/(\\d{4}).", "\\1", taxa_count)),
                          grepl("^\\d{4}[:]", taxa_count) ~ sub("[:].*", "", taxa_count),
                          T ~ date),
         taxa_count = case_when(taxa_count == "NULL" ~ NA,
                                grepl("no data", taxa_count, ignore.case = TRUE) ~ NA,
                                grepl("^in ", taxa_count, ignore.case = TRUE) ~ sub(".*[;,] ", "", taxa_count),
                                grepl("information not provided", taxa_count, ignore.case = TRUE) ~ NA,
                                grepl("date to year only", taxa_count, ignore.case = TRUE) ~ sub("(\\d+) observed[.;].*", "\\1", taxa_count),
                                grepl("number of specimens", taxa_count, ignore.case = TRUE) ~ sub("Number of Specimens: ", "", taxa_count),
                                grepl("one fish collected|one fish present", taxa_count, ignore.case = TRUE) ~ "1",
                                grepl("unpublished data", taxa_count, ignore.case = TRUE) ~ gsub("[^0-9]", "", taxa_count),
                                grepl("early season, positive|late season, positive", taxa_count, ignore.case = TRUE) ~ "1",
                                grepl("^\\d{4}:", taxa_count) ~ sub("^.*:", "", taxa_count),
                                T ~ taxa_count),
         taxa_count = case_when(grepl("observed$|observerd$|present$|Adult$|collected$|individuals$", taxa_count, ignore.case = TRUE) ~ sub("[A-Za-z\\s]+", "", taxa_count),
                                grepl("creek|river|lake|pits|fork", taxa_count, ignore.case = TRUE) ~ sub("[,.].*", "", taxa_count),
                                T ~ taxa_count),
         taxa_count = case_when(grepl("^-", taxa_count) ~ sub(".*\\b(\\d+)\\b(?=\\s+present).*", "\\1", taxa_count, perl = TRUE),
                                grepl("WERE COLLECTED", taxa_count) ~ "1",
                                T ~ taxa_count),
         taxa_count = sub("^\\s", "", taxa_count),
         date = case_when(grepl("^\\d{4}[;]", taxa_count) ~ sub("[;].*", "", taxa_count),
                          T ~ date),
         taxa_count = sub(".*\\b(\\d+)\\s*(specimen|individual|collected|adult|observed|male|female|paratype|darters observed|indiv|young|seen|sighted|breed|dead|leopard|larvae|ad|ind).*", "\\1", taxa_count, ignore.case = T),
         taxa_count = case_when(grepl("one", taxa_count, ignore.case = TRUE) ~ "1",
                                grepl("two", taxa_count, ignore.case = TRUE) ~ "2",
                                grepl("three", taxa_count, ignore.case = TRUE) ~ "3",
                                grepl("four", taxa_count, ignore.case = TRUE) ~ "4",
                                grepl("five", taxa_count, ignore.case = TRUE) ~ "5",
                                grepl("six", taxa_count, ignore.case = TRUE) ~ "6",
                                grepl("seven", taxa_count, ignore.case = TRUE) ~ "7",
                                grepl("nine", taxa_count, ignore.case = TRUE) ~ "9",
                                grepl("eight", taxa_count, ignore.case = TRUE) ~ "8",
                                grepl("twelve", taxa_count, ignore.case = TRUE) ~ "12",
                                grepl("twenty-six", taxa_count, ignore.case = TRUE) ~ "26",
                                grepl("no ", taxa_count, ignore.case = TRUE) ~ "0",
                                grepl("^\\d{4};", taxa_count) ~ sub("present", "", sub("^.*; ", "", taxa_count)),
                                grepl("female found|male found|found while|and released|pigg|mm|duke|UV|black|river|odwc|creek", taxa_count, ignore.case = TRUE) ~ "1",
                                T ~ taxa_count),
         taxa_count = ifelse(grepl("^[^0-9]*$", taxa_count), "1", taxa_count),
         taxa_count = ifelse(grepl("spec", taxa_count, ignore.case = TRUE), "6", taxa_count),
         taxa_count = ifelse(grepl("captured and relocated", taxa_count), "2", taxa_count),
         taxa_count = ifelse(grepl("^Crew", taxa_count), "20", taxa_count),
         taxa_count = ifelse(taxa_count == "1983", "1", taxa_count),
         taxa_count = ifelse(grepl("and", taxa_count), regmatches(taxa_count, regexpr("\\d+", taxa_count)), taxa_count),
         taxa_count = as.numeric(taxa_count)) %>%
  group_by(lat, long) %>%
  mutate(site_id = paste("ONHI-", cur_group_id(), sep = "")) %>%
  ungroup() %>%
  distinct()

write_csv(onhi, paste0(PATH, "/01_BioDat/source_dat_cleaned/onhi_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(onhi)


#QUINN ----
## convert shp to csv
shp.list <- list.files(dat_path, pattern = "Quinn.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

lapply(shp.list, function(x) {
  st_read(dsn = dat_path, layer = x) %>%
    st_drop_geometry() %>%
    write_csv(., file = paste0(dat_path, "/", x, ".csv"))
})

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "Quinn.*\\.csv$", full.names = T) %>% unique()
lapply(file.list, function(x) head(read_csv(x)))

quinn <- list()
for(i in 2:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  quinn[[i]] <- fix_clean_occ(file.list[[i]],
                              species.id.col = "species", 
                              site.id.col = "site",
                              set.source = "Quinn_AGFC", save.it = NA)
  
  flush.console()
}

q1 <- read_csv(file.list[[1]]) %>%
  rename(lat = Lat_dd, long = Lon_dd, site_id = 'Site') %>%
  mutate(species = "Percina pantherina",
         across(.cols = contains("19")|contains("20"), ~ as.numeric(gsub("\\D", "", .x)))) %>%
  pivot_longer(cols = contains("19")|contains("20"), names_to = "date", values_to = "taxa_count") %>%
  select(site_id, lat, long, species, date, taxa_count) %>%
  filter(!is.na(taxa_count)) %>%
  distinct() %>%
  mutate(bio_type = "fish", source = "Quinn_AGFC")

quinn <- bind_rows(quinn) %>% distinct() %>% mutate(bio_type = "fish")

quinn <- bind_rows(quinn, q1) %>% distinct()

quinn <- quinn %>%
  group_by(across(-c(date, taxa_count))) %>%
  filter(n() > 1, !is.na(date)) %>% #remove rows that are NA in the date column, if that group already has other data
  ungroup()

write_csv(quinn, paste0(PATH, "/01_BioDat/source_dat_cleaned/quinn_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(quinn, q1)

#RAM ----
readxl::read_xlsx(paste0(dat_path, "/RAM_Invert_FOX_UTMs.xlsx")) %>% write_csv(., paste0(dat_path, "/RAM_Invert_FOX_UTMs_2.csv"))
file.list <- list.files(dat_path, pattern = "RAM.*\\.csv$", full.names = T) %>% unique()
lapply(file.list, function(x) head(read_csv(x)))

ram <- lapply(file.list, function(x) {
  df <- read_csv(x) %>%
    select(Site, UTM_X, UTM_Y, Taxa, Totals, Date_col) %>%
    rename(site_id = Site, taxa_name = Taxa, taxa_count = Totals, date = Date_col) %>%
    mutate(date = as.character(date))
  return(df)
})

ram <- lapply(ram, function(x){
  err <- try(mutate(x, date = as.Date(date)))
  if(inherits(err, "try-error")) {
    x <- mutate(x, date = as.Date(format(strptime(date, format = "%d-%b-%y"), "%Y-%m-%d")))
  } else {
    x <- mutate(x, date = as.Date(date))
  }
})

ram <- bind_rows(ram) %>% distinct() %>% filter(!is.na(UTM_X) & !is.na(UTM_Y))
r2 <- ram %>%
  st_as_sf(., coords = c("UTM_X", "UTM_Y"), crs = 26915, remove = FALSE) %>%
  st_transform(., crs = 4326)
ram <- cbind(ram, st_coordinates(r2)) %>% 
  select(-UTM_X, -UTM_Y) %>%
  rename(lat = Y, long = X) %>%
  mutate(bio_type = "bug", source = "RAM")
ram <- ram %>% mutate(lat = ifelse(lat < 10, lat * 10, lat)) #fix column coordinates

write_csv(ram, paste0(PATH, "/01_BioDat/source_dat_cleaned/ram_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(ram, r2)

#REQUEST UNKNOWN SRCE ----
file.list <- list.files(dat_path, pattern = "Request.*\\.xlsx$", full.names = T) %>% unique()
lapply(file.list, function(x){
  readxl::read_xlsx(x) %>% 
    write_csv(., paste0(sub(".xlsx", "", x), ".csv"))
})

file.list <- list.files(dat_path, pattern = "Request.*\\.csv$", full.names = T) %>% unique()
lapply(file.list, function(x) head(read_csv(x)))

ok_unk <- list()
for(i in 1:length(file.list)) {
  cat("\n\nCleaning", file.list[[i]], "\n")
  ok_unk[[i]] <- fix_clean_occ(file.list[[i]], ignore.parse.err = TRUE,
                               site.id.col = "WBID",
                               count.id.col = c("Total", "NumRiff", "NumVeg", "NumWoody"),
                               set.source = "OK_UNK", save.it = NA)
  
  flush.console()
}

ok_unk <- bind_rows(ok_unk)

ok_unk <- ok_unk %>%
  filter(!if_all(c("order", "family", "genus", "species"), is.na)) %>%
  filter(!if_any(c("order", "family", "genus", "species"), ~grepl("hybrid", .x))) %>%
  distinct()

write_csv(ok_unk, paste0(PATH, "/01_BioDat/source_dat_cleaned/OK_UNK_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(ok_unk)


#SMB ----
shp.list <- list.files(dat_path, pattern = "SMB.*\\.shp", full.names = F) %>%
  sub("\\.shp.*$", "", .) %>%
  unique()

lapply(shp.list, function(x) {
  st_read(dsn = dat_path, layer = x) %>%
    st_drop_geometry() %>%
    write_csv(., file = paste0(dat_path, "/", x, ".csv"))
})

#clean data + bind + save
file.list <- list.files(dat_path, pattern = "SMB.*\\.csv$", full.names = T) %>% unique()
SMB <- lapply(file.list, function(x) {
  df <- read_csv(x) %>%
    rename(site_id = ReachCode, taxa_count = M_dolomieu) %>%
    mutate(species = "Micropterus dolomieu")
  return(df)
  })

SMB <- bind_rows(SMB) %>% distinct()
r2 <- SMB %>%
  st_as_sf(., coords = c("site_x", "site_y"), crs = 26915, remove = FALSE) %>%
  st_transform(., crs = 4326)
SMB <- cbind(SMB, st_coordinates(r2)) %>% 
  select(-site_x, -site_y) %>%
  rename(lat = Y, long = X) %>%
  mutate(bio_type = "fish", source = "SMB")

SMB <- SMB %>% 
  group_by(lat, long) %>%
  mutate(site_id = paste("SMB-", cur_group_id(), sep = "")) %>%
  ungroup() %>%
  distinct()

write_csv(SMB, paste0(PATH, "/01_BioDat/source_dat_cleaned/SMB_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(SMB, r2)

#University of AR ----
readxl::read_xlsx(paste0(dat_path, "/University of AR Tyler Fox fish and bugs.xlsx"),
                         sheet = 2) %>%
  write_csv(., paste0(dat_path, "/UnivAR_TylerFox_bugs.csv"))
file.list <- list.files(dat_path, pattern = "UnivAR.*\\.csv$", full.names = T) %>% unique()

ua_bug <- read_csv(paste0(dat_path, "/UnivAR_TylerFox_bugs.csv")) %>%
  select('Activity Start Date', 'Monitoring Location ID', 'Monitoring Location Latitude', 'Monitoring Location Longitude',
         'Characteristic Name', 'Result Value', 'Result Unit', 'Taxonomic Name') %>%
  rename(date = 'Activity Start Date', site_id = 'Monitoring Location ID', lat = 'Monitoring Location Latitude',
         long = 'Monitoring Location Longitude', taxa_name = 'Taxonomic Name', taxa_count = 'Result Value') %>%
  filter(`Characteristic Name` == "Count") %>%
  select(-`Characteristic Name`, -`Result Unit`) %>%
  distinct() %>%
  mutate(taxa_count = ifelse(taxa_count == 0, NA, taxa_count)) %>%
  group_by(across(-c(date, taxa_count))) %>%
  filter(n() > 1, !is.na(taxa_count)) %>% #remove rows that are NA in the date column, if that group already has other data
  ungroup() %>%
  mutate(bio_type = "bug", source = "OWSB",
         taxa_name = ifelse(grepl("retire", taxa_name, ignore.case = TRUE), sub(".*use ", "", taxa_name), taxa_name))

ua_fish <- read_csv(paste0(dat_path, "/UnivAR_TylerFox_fish.csv")) %>%
  select('Activity Start Date', 'Monitoring Location ID', 'Monitoring Location Latitude', 'Monitoring Location Longitude',
         'Characteristic Name', 'Result Value', 'Taxonomic Name') %>%
  rename(date = 'Activity Start Date', site_id = 'Monitoring Location ID', lat = 'Monitoring Location Latitude',
         long = 'Monitoring Location Longitude', species = 'Taxonomic Name', taxa_count = 'Result Value') %>%
  filter(`Characteristic Name` == "Count") %>%
  select(-`Characteristic Name`) %>%
  distinct() %>%
  mutate(bio_type = "fish", source = "OWSB",
         species = ifelse(grepl("retire", species, ignore.case = TRUE), sub(".*use ", "", species), species))

write_csv(bind_rows(ua_bug, ua_fish), paste0(PATH, "/01_BioDat/source_dat_cleaned/OWSB2_alltax_clean_", gsub("-", "", Sys.Date()), ".csv"))
rm(ua_bug, ua_fish)

occ.dat <- lapply(file.list, function(x) read_csv(x, col_types = cols(site_id = col_character(), date = col_character(), date_min = col_character(), date_max = col_character())))

#compiled data, final fixes ----
occ.dat <- bind_rows(occ.dat) %>% distinct()

occ.dat <- occ.dat %>%
  mutate(lat2 = ifelse(lat < 0, long, lat),
         long2 = ifelse(lat < 0, lat, long),
         lat2 = ifelse(lat2 < 10, lat2 * 10, lat2)) %>%
  select(-lat, -long, -HybridFlag) %>% #already been removed
  rename(lat = lat2, long = long2) %>%
  mutate(long = abs(long)*-1,
         source = ifelse(is.na(source), "Matthews", source)) %>%
  relocate(site_id, lat, long, date, order, family, genus, species, taxa_name, date_min, date_max, date_equal) %>%
  distinct()

write_csv(occ.dat, paste0(PATH, "/01_BioDat/occ_alltaxa_combined_nofilt_", gsub("-", "", Sys.Date()), ".csv"))




