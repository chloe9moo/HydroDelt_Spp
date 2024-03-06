#### FUNCTIONS FOR TRAIT PREP ####

#for finding most common categorical in na replacement
find_cat_mode <- function(x) { 
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#for replacing NAs in traits with the mean for contin. or mode for categorical
replace_na_traits <- function(data, group_var, cat_vars) {
  data %>%
    { if (!is.na(group_var)) group_by(., !!sym(group_var)) else . } %>%
    mutate(across(c(where(is.numeric), -matches("record_ttl_15kbuff")), 
                  ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
    mutate(across(matches(paste0(cat_vars, collapse = "|")), 
                  ~ ifelse(is.na(.x), find_cat_mode(.x[!is.na(.x)]), .x))) %>%
    ungroup()
}

#for summarizing to one value per species
compress_traits <- function(data, group_var, cat_vars) {
  data %>%
    { if (!is.na(group_var)) group_by(., !!sym(group_var)) else . } %>%
    summarise(across(c(where(is.numeric), -matches("record_ttl_15kbuff")), 
                     ~ mean(.x, na.rm = TRUE)),
              across(matches(paste0(cat_vars, collapse = "|")), 
                     ~ find_cat_mode(.x[!is.na(.x)]))) %>%
    ungroup()
}

#make site x trait matrix
site_x_trait <- function(site_dat, #wide format site data
                         trait_dat, #trait data, no NAs, only species name and traits
                         tax_col_name = c("taxa", "species"), #column location for taxa ID name
                         site_cols = c("lat", "long", "COMID", "flw_type", "gage_no_15yr", "dist2gage_m_15yr", "dist2strm_m_flw", "site_id"), #site cols to ignore
                         convert_cont_to_cat = TRUE,
                         matrix_type = c("presence/absence", "abundance")) 
{
  #prep site matrix for joining to traits
  o <- site_dat %>%
    mutate(tot_tax = rowSums(select(., -all_of(site_cols)))) %>% #for later
    pivot_longer(cols = -c(tot_tax, any_of(site_cols)), names_to = "taxa", values_to = "p") %>%
    filter(p == 1) %>%
    select(-p)
  
  #prep traits for joining to site matrix
  t <- trait_dat %>%
    mutate(across(.cols = any_of(tax_col_name), ~ gsub(" ", "_", .x)),
           across(.cols = where(is.character), ~ gsub("/", "-", .x)))
  
  ##get already binary columns (that aren't categorical), and continuous vars
  bi_col <- sapply(t, function(x) length(unique(na.omit(x))) == 2 & !is.character(x))
  bi_col <- names(bi_col[bi_col == TRUE])
  
  ##get empty col
  zero_col <- sapply(t, function(x) ifelse(is.numeric(x), ifelse(sum(x) == 0, TRUE, FALSE), FALSE))
  zero_col <- names(zero_col[zero_col == TRUE])
  t <- t %>% select(-all_of(zero_col))
  
  ##categorize continuous variables
  if(convert_cont_to_cat) {
    t <- t %>%
      mutate(across(.cols = c(where(is.numeric), -all_of(bi_col)),
                    ~ case_when(.x <= quantile(.x, c(0.25)) ~ "low",
                                .x > quantile(.x, c(0.25)) & .x <= quantile(.x, c(0.50)) ~ "low-med",
                                .x > quantile(.x, c(0.50)) & .x <= quantile(.x, c(0.75)) ~ "med-high",
                                .x > quantile(.x, c(0.75)) ~ "high")))
  }
  
  ##make traits long format for joining
  t <- t %>%
    pivot_longer(cols = -c(any_of(tax_col_name), all_of(bi_col)), names_to = "trait", values_to = "value") %>%
    mutate(p = 1) %>%
    pivot_wider(names_from = c(trait, value), values_from = p, values_fill = list(p = 0)) %>%
    pivot_longer(cols = -any_of(tax_col_name), names_to = "trait", values_to = "p") %>%
    filter(p == 1)
  
  #join site x traits
  if(!"taxa" %in% names(t)) { #rename taxa name column for joining
    names(t)[names(t) %in% tax_col_name] <- "taxa"
  }
  
  #check for matching taxa names
  if(!all(o$taxa %in% t$taxa)) {
    stop("Trait taxa and site taxa do not match!")
  }
  
  t.o <- left_join(o, t, relationship = "many-to-many") %>%
    select(-any_of(tax_col_name)) %>%
    pivot_wider(names_from = trait, values_from = p, values_fn = ~ sum(.x, na.rm = FALSE))
  
  if(nrow(site_dat) != nrow(t.o)) { message("Oops! Something went wrong. Original # of sites and # of sites in new matrix are not the same.") }
  if(!any(c("presence/absence", "abundance") %in% matrix_type)) { stop("Unsupported matrix type specified. Check spelling.") }
  
  out_list <- list()
  
  if("presence/absence" %in% matrix_type) { #get presence/absence matrix (any species with modality = present)
    pa <- t.o %>%
      select(-tot_tax) %>%
      mutate(across(-any_of(site_cols),
                    ~ ifelse(is.na(.x), 0, 1)))
    out_list[[1]] <- pa 
    names(out_list)[1] <- "presence_absence"
  }
  
  if("abundance" %in% matrix_type) {
    abund <- t.o %>%
      mutate(across(.cols = -c(tot_tax, any_of(site_cols)), 
                    ~ ifelse(is.na(.x), 0, .x / tot_tax))) %>%
      select(-tot_tax)
    out_list[[2]] <- abund
    names(out_list)[2] <- "abundance"
  }
  
  #remove empty list elements if applicable
  indx <- !sapply(out_list, is.null)
  out_list <- out_list[indx]
  
  #return adjusted matrix(ces)
  if(length(out_list) == 1) {
    return(out_list[[1]])
  } else {
    return(out_list)
  }
  
}
