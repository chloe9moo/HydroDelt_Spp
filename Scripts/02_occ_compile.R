## Compile and fix species names occurrence data for fish and inverts

library(tidyverse); library(rgbif)#; library(sf)
options(readr.show_col_types = FALSE)

PATH <- getwd()

dat_path <- paste0(PATH, "/01_BioDat/source_dat_cleaned")

file.list <- list.files(dat_path, ".csv", full.names = TRUE)

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

#get and fix all taxa names ----
spp_names <- occ.dat %>%
  select(order, family, genus, species, taxa_name) %>%
  distinct() %>%
  mutate(key = row_number())
search_nam <- spp_names %>%
  pivot_longer(-key, names_to = "level", values_to = "name") %>%
  distinct() %>%
  filter(!is.na(name))

# crul::set_opts(http_version = 2)

k_ranks <- search_nam %>% filter(level != "taxa_name") %>%
  group_by(level, name) %>%
  summarise(key = toString(key)) %>%
  ungroup()
u_ranks <- search_nam %>% filter(level == "taxa_name") %>%
  group_by(level, name) %>%
  summarise(key = toString(key)) %>%
  ungroup()

taxa_nam <- NULL
for(i in 1:nrow(k_ranks)) {
  taxa1 <- name_backbone(name = k_ranks[i, "name"], rank = k_ranks[i, "level"], kingdom = "animals")
  taxa_nam <- bind_rows(taxa_nam, taxa1)
}

k_ranks <- cbind(k_ranks, taxa_nam)

missing_nam <- k_ranks %>% filter(matchType == "NONE") %>%
  select(level, name, key) %>%
  mutate(new_name = case_when(name == "Crooked Creek" ~ "remove",
                              name == "Acariformes" ~ "Acariformes",
                              name == "Ceriantharia" ~ "Ceriantharia",
                              name == "Cladocera" ~ "Cladoceromorpha",
                              name == "Diamesinae" ~ "Diamesinae",
                              name == "Ermyzon" ~ "Erimyzon",
                              name == "Exos niger" ~ "Esox niger",
                              name == "Hydroida" ~ "Hydroidolina",
                              name == "Littoridinomorpha" ~ "Littorinimorpha",
                              name == "Lumbricina" ~ "Crassiclitellata",
                              name == "Mesogastropoda" ~ "Caenogastropoda",
                              name == "Neotaenioglossa" ~ "Littorinida",
                              name == "Orthocladiinae" ~ "Orthocladiinae",
                              name == "Paleoheterodonta" ~ "Heteroconchia",
                              name == "Pharyngobdellida" ~ "Arhynchobdellida",
                              name == "Pseudochironomini" ~ "Pseudochironomus", #pseudochrironomini is a tribe with only one genus, so combining
                              name == "Scolecida" ~ "Scolecida",
                              name == "Semionotiformes" ~ "Lepisosteiformes", #the original is an extinct genus (fossils)
                              name == "Tanypodinae" ~ name,
                              name == "Chironomini" ~ "Chironominae", #this is a tribe, collapsing into the subfamily
                              name == "Tanytarsini" ~ "Chironominae", #this is a tribe, collapsing into the subfamily
                              name == "Unidentified" ~ "remove",
                              name == "Unidentified spp." ~ "remove",
                              name == "Unidentified YOY" ~ "remove",
                              name == "Unidentified Etheostoma spp." ~ "Etheostoma",
                              name == "Unidentified Fundulus spp." ~ "Fundulus",
                              name == "Unidentified Ictiobus spp." ~ "Ictiobus",
                              name == "Unidentified Lepisosteus spp." ~ "Lepisosteus",
                              name == "Unidentified Lepomis spp." ~ "Lepomis",
                              name == "Unidentified Notropis spp." ~ "Notropis",
                              name == "Unidentified Percidae spp." ~ "Percidae",
                              name == "Unidentified Phoxinus spp." ~ "Phoxinus",
                              name == "Unidentified Micropterus spp." ~ "Micropterus",
                              name == "Unidentified Pimephales spp." ~ "Pimephales",
                              name == "Unionacea" ~ "Unionoidea",
                              grepl("bullhead", name) ~ "Ameiurus",
                              grepl("Catostomidae", name) ~ "Catostomidae",
                              grepl("hybrid", name) ~ "remove",
                              grepl(" x ", name) ~ "remove",
                              grepl("acari|arcari", name, ignore.case = TRUE) & name != "Acariformes" ~ "Arachnida", #in gbif backbone, match to that
                              grepl("ammoc", name, ignore.case = TRUE) ~ "Petromyzontiformes",
                              grepl("attract", name, ignore.case = TRUE) ~ "Atractosteus",
                              grepl("Basommatophora", name) ~ "Heterobranchia",
                              grepl("lythrurus", name, ignore.case = TRUE) ~ "Lythrurus",
                              grepl("moxostoma", name, ignore.case = TRUE) ~ "Moxostoma"))
df <- data.frame(name = c("Acariformes", "Ceriantharia", "Cladocera", "Diamesinae", "Hydroida", "Mesogastropoda", "Neotaenioglossa", "Orthocladiinae", "Paleoheterodonta", "Scolecida", "Tanypodinae", "Chironomini", "Tanytarsini", "Unionacea"), 
           new_name = c("Acariformes", "Ceriantharia", "Cladoceromorpha", "Diamesinae", "Hydroidolina", "Caenogastropoda", "Littorinida", "Orthocladiinae", "Heteroconchia", "Scolecida", "Tanypodinae", "Chironominae", "Chironominae", "Unionoidea"),
           level = c("superorder", "order", "infraorder", "subfamily", "subclass", "infraclass", "order", "subfamily", "infraclass", "infraclass", "subfamily", "subfamily", "subfamily", "superfamily"),
           kingdom = rep("Animalia", 14),
           phylum = c("Arthropoda", "Cnidaria", "Arthropoda", "Arthropoda", "Cnidaria", "Mollusca", "Mollusca", "Arthropoda", "Mollusca", "Annelida", "Arthropoda", "Arthropoda", "Arthropoda", "Mollusca"),
           class = c("Arachnida", "Anthozoa", "Branchiopoda", "Insecta", "Hydrozoa", "Gastropoda", "Gastropoda", "Insecta", "Bivalvia", "Polychaeta", "Insecta", "Insecta", "Insecta", "Bivalvia"), 
           order = c(NA, "Ceriantharia", "Diplostraca", "Diptera", NA, NA, NA, "Diptera", NA, NA, "Diptera", "Diptera", "Diptera", "Unionida"),
           family = c(NA, NA, NA, "Chironomidae", NA, NA, NA, "Chironomidae", NA, NA, "Chironomidae", "Chironomidae", "Chironomidae", NA))

missing_nam <- left_join(missing_nam, df, by = c("new_name", "name"))



#get all the unique codes to search for u_ranks
taxa_key <- k_ranks %>% 
  select(kingdom, phylum, class, order, family, genus, species, 
         kingdomKey, phylumKey, classKey, orderKey, familyKey, genusKey, speciesKey) %>%
  distinct() %>%
  pivot_longer(., where(is.numeric), names_to = "name_key", values_to = "name_key_value") %>%
  mutate(name_key = sub("Key", "", name_key),
         name_key_result = case_when(name_key == "kingdom" ~ kingdom,
                                     name_key == "phylum" ~ phylum,
                                     name_key == "class" ~ class,
                                     name_key == "order" ~ order,
                                     name_key == "family" ~ family,
                                     name_key == "genus" ~ genus,
                                     name_key == "species" ~ species)) %>%
  select(contains("name_key")) %>% distinct() %>% filter(!is.na(name_key_result))

u_ranks %>%
  left_join(., taxa_key, by = c("name" = "name_key_result"))

missing_nam %>%
  left_join(., taxa_key, by = c("new_name" = "name_key_result")) %>% View()
    
k_ranks %>%
  separate_longer_delim(., key, delim = ",") %>% View()


#remove synonyms with repeat species keys
pnw.gbif.taxa <- pnw.gbif.taxa %>% filter(!duplicated(speciesKey))




