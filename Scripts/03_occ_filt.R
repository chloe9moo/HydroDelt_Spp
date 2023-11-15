# OCCURRENCE FILTERING

library(tidyverse)

PATH <- getwd()

occ.dat <- read_csv(paste0(PATH, "/01_BioDat/occ_alltax_updatedspp_nofilt_20231113.csv"),  col_types = cols(date_min = col_character(), date_max = col_character(),
                                                                                                            subclass = col_character(), superorder = col_character(), 
                                                                                                            suborder = col_character(), infraorder = col_character(), 
                                                                                                            superfamily = col_character(), subgenus = col_character(), 
                                                                                                            species = col_character(), subspecies = col_character()))

occ.dat <- occ.dat %>%
  select(bio_site_no, lat, long, contains("date"), taxa_count,
         kingdom, phylum, 
         class, subclass,
         superorder, order, suborder, infraorder,
         superfamily, family, subfamily, 
         tribe, genus, subgenus, species, subspecies,
         -contains("orig")) %>%
  mutate(taxa_count = ifelse(is.na(taxa_count), 1, taxa_count)) %>%
  filter(kingdom != "remove") %>% #get rid of those flagged for removal
  filter(taxa_count != 0) #get rid of zero count ids

table(occ.dat$phylum)

#get fish dat
fish <- occ.dat %>%
  filter(phylum == "Chordata") %>% #verts
  filter(!class %in% c("Aves", "Amphibia") & order != "Testudines") #remove non-fish verts

write_csv(fish, paste0(PATH, "/01_BioDat/occ_fish_all_long_20231113.csv"))

#get invert dat
bug <- occ.dat %>%
  filter(class == "Insecta")

write_csv(fish, paste0(PATH, "/01_BioDat/occ_bug_all_long_20231113.csv"))
