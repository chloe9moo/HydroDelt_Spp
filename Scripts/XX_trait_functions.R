#### FUNCTIONS FOR TRAIT PREP ####

# #for finding most common categorical in na replacement
find_cat_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# #for replacing NAs in traits with the mean for contin. or mode for categorical
# replace_na_traits <- function(data, group_var, cat_vars) {
#   data %>%
#     { if (!is.na(group_var)) group_by(., !!sym(group_var)) else . } %>%
#     mutate(across(c(where(is.numeric), -matches("record_ttl_15kbuff")), 
#                   ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
#     mutate(across(matches(paste0(cat_vars, collapse = "|")), 
#                   ~ ifelse(is.na(.x), find_cat_mode(.x[!is.na(.x)]), .x))) %>%
#     ungroup()
# }

#for summarizing to one value per species
# compress_traits <- function(data, group_var, cat_vars) {
#   data %>%
#     { if (!is.na(group_var)) group_by(., !!sym(group_var)) else . } %>%
#     summarise(across(c(where(is.numeric), -matches("record_ttl_15kbuff")), 
#                      ~ mean(.x, na.rm = TRUE)),
#               across(matches(paste0(cat_vars, collapse = "|")), 
#                      ~ find_cat_mode(.x[!is.na(.x)]))) %>%
#     ungroup()
# }

#GET TRAIT GROUP CLUSTERS
#FOR TESTING
# site_dat <- occ.list[[1]]
# trait_dat <- traits[[1]][,-c(1:3)]

# library(cluster); library(vegan)
#prep trait data for distance matrix
prep_trait_df_cols <- function(trait_dat,
                               tax_col_name = "species",
                               ord_col # list containing ordinal traits, in the order specified (from first to last)
                               ) {
  #set names
  trait_dat <- trait_dat %>%
    { if("species" %in% names(.)) rename(., taxa = species) else . } %>%
    mutate(taxa = gsub(" ", "_", taxa)) %>%
    column_to_rownames("taxa")
  
  #set up binary, ordinal, categorical traits as factors as appropriate
  bi_col <- sapply(trait_dat, function(x) length(unique(na.omit(x))) == 2 & !is.character(x))
  bi_col <- names(bi_col[bi_col])
  
  cat_col <- sapply(trait_dat, function(x) is.character(x))
  cat_col <- names(cat_col[cat_col])
  cat_col <- cat_col[!cat_col %in% names(ord_col)]
  
  trait_dat <- trait_dat %>%
    mutate(across(matches(c(bi_col, cat_col)), ~ as.factor(.x)))
  
  for (col_name in names(ord_col)) {
    if (col_name %in% names(trait_dat)) {
      trait_dat[[col_name]] <- factor(trait_dat[[col_name]], levels = ord_col[[col_name]], ordered = TRUE)
    }
  }
  
  return(trait_dat)
}

#silhouette method for finding most likely num of clusters
id_best_grp <- function(k, h_clust, dist_mat, return_clust = FALSE) { #get avg. silhouette for k clusters
  message("Testing ", k, " groups..")
  
  grps <- cutree(h_clust, k = k)
  ss <- silhouette(grps, dist_mat)
  mn.ss <- mean(ss[, 3])
  
  #significance of groups with PERMANOVA
  grps <- factor(grps, levels = seq(1, k))
  ado <- adonis2(dist_mat ~ grps, permutations = 999)
  
  #get results
  tmp <- vector("list", 0)
  tmp[["mean_silhouette_width"]] <- mn.ss
  tmp[["permanova"]] <- ado
  
  if(return_clust) {
    tmp[["cluster_assignment"]] <- grps
  }
  
  return(tmp)
}

comp_trait_groups <- function(trait_dat, #df with species name column, and then all other traits
                              hc_method = c("average", "single", "complete", "ward"), #one or try multiple
                              max_k = 100, #maximum number of groups to test and compare,
                              tax_col_name = "species",
                              ord_col #set other variables for internal functions
                              ) {
  #prep trait df
  trait_dat <- prep_trait_df_cols(trait_dat = trait_dat, ord_col = ord_col)
  
  #gower distance matrix
  g.dist <- daisy(trait_dat, metric = "gower")
  
  #hierarchical clustering
  # dhc <- diana(g.dist)
  if(length(hc_method) > 1) { #compare multiple clustering methods and find highest agglomerative coeff
    names(hc_method) <- hc_method
    
    tmp <- sapply(hc_method, function(x) agnes(g.dist, method = x)$ac) #function to get coefficient
    
    hc_method <- hc_method[which.max(tmp)]
  }

  ahc <- agnes(g.dist, method = hc_method)
  
  #get groups
  # Compute and plot silhouette width for k=2 to number of traits (don't think it makes sense to have more than that)
  clust.n <- 2:max_k
  
  # extract avg silhouette for 2-15 clusters
  clust.test <- lapply(clust.n, id_best_grp, h_clust = ahc, dist_mat = g.dist)
  
  clust.df <- data.frame()
  for(i in seq_along(clust.test)) {
    tmp <- data.frame(k = i, 
               mn_sil_w = clust.test[[i]]$mean_silhouette_width,
               grp_r2 = clust.test[[i]]$permanova[1, "R2"],
               F_stat = clust.test[[i]]$permanova[1, "F"], 
               p_val = clust.test[[i]]$permanova[1, "Pr(>F)"])
    clust.df <- bind_rows(clust.df, tmp)
  }
  
  output_list <- list()
  
  output_list[["hierarchical_cluster_results"]] <- ahc
  output_list[["gower_matrix"]] <- g.dist
  output_list[["group_comp_results"]] <- clust.df
  
  return(output_list)
}

#pcoa for looking at various cluster numbers
plot_clusters <- function(
    trait_group_output, #results returned from 'comp_trait_groups' function
    num_clust, #number of clusters to examine
    return_pc, #vector of length 2, for which PCs to plot
    ellipse = TRUE, #include group ellipses in plot?
    vectors = TRUE, #include trait vectors in plot?
    clust_legend = FALSE, #include legend for cluster colors?
    #only necessary if plotting vectors (needed to get them)
    trait_dat = traits[[1]][,-c(1:3)],
    ord_col = oc
) {
  
  pcoa <- cmdscale(trait_group_output[["gower_matrix"]], k = num_clust, eig=TRUE)
  
  var_exp <- (pcoa$eig/sum(pcoa$eig)*100)[return_pc] #variation explained by each PC
  var_exp <- round(var_exp, digits = 2)
  
  #get species groups
  grps <- cutree(trait_group_output[["hierarchical_cluster_results"]], k = num_clust) #get groups
  
  #get species coordinates + group numbers
  suppressMessages(
    pc.scores <- bind_cols(row.names(pcoa$points), pcoa$points[, return_pc], grps) #add scores to traits dataset
  )
  names(pc.scores) <- c("taxa", paste0("PC", return_pc), "cluster")
  pc.scores$cluster <- factor(pc.scores$cluster, levels = seq(1, num_clust))
  
  #get trait arrows
  if(vectors) {
    tr.df <- prep_trait_df_cols(trait_dat = trait_dat, ord_col = ord_col)
    vec.trait <- envfit(scores(pcoa), tr.df, perm=1000, na.rm = TRUE)
    t.scrs <- as.data.frame(scores(vec.trait, display = "vectors"))
    t.scrs <- cbind(t.scrs, trait = rownames(t.scrs))
  }
  
  p1 <- ggplot() +
    { if(ellipse) { stat_ellipse(data = pc.scores, geom = "polygon", 
                                 aes(.data[[names(pc.scores)[2]]], .data[[names(pc.scores)[3]]], color = cluster), 
                                 alpha = 0.5, fill = NA, linewidth = 2, show.legend = F, type = "norm", level = 0.9)}} +
    geom_point(data = pc.scores, aes(.data[[names(pc.scores)[2]]], .data[[names(pc.scores)[3]]], fill = cluster), shape = 21, show.legend = clust_legend) +
    { if(vectors) { geom_segment(data = t.scrs,
                                 aes(x = 0, xend = Dim1/2, y = 0, yend = Dim2/2),
                                 arrow = arrow(length = unit(0.25, "cm")), colour = "black")}} +
    { if(vectors) { geom_text(data = t.scrs, aes(x = Dim1/1.5, y = Dim2/1.5, label = trait), size = 3, angle = 0)}} +
    xlab(paste0(names(pc.scores)[2]," (", round(var_exp[1], digits = 2), " %)" )) +
    ylab(paste0(names(pc.scores)[3]," (", round(var_exp[2], digits = 2), " %)" )) +
    theme_minimal()
  
  return(p1)
}

# #make site x trait matrix
site_x_trait <- function(site_dat, #wide format site data
                         trait_dat, #trait data, no NAs, only species name and traits
                         tax_col_name = c("taxa", "species"), #column location for taxa ID name
                         ignore_cols = c("lat", "long", "COMID", "flw_type", 
                                       "gage_no_15yr", "dist2gage_m_15yr", "dist2strm_m_flw", "site_id"), #site cols to ignore
                         convert_cont_to_cat = TRUE, #convert continuous traits to categorical
                         match_spp_across_dat = TRUE, #filter spp to match within each dataset (spp in traits will = spp in occ, LCD type of thing)
                         matrix_type = c("presence/absence", "abundance")) 
{
  #prep site matrix for joining to traits
  o <- site_dat %>%
    mutate(tot_tax = rowSums(select(., -all_of(ignore_cols)))) %>% #for later
    pivot_longer(cols = -c(tot_tax, any_of(ignore_cols)), names_to = "taxa", values_to = "p") %>%
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
  
  ##get continuous cols
  cont_col <- select(t, c(where(is.numeric), -all_of(bi_col)))
  cont_col <- names(cont_col)
  
  ##categorize continuous variables
  if(convert_cont_to_cat) {
    t <- t %>%
      mutate(across(.cols = matches(cont_col),
                    ~ case_when(.x <= quantile(.x, c(0.25)) ~ "low",
                                .x > quantile(.x, c(0.25)) & .x <= quantile(.x, c(0.50)) ~ "low-med",
                                .x > quantile(.x, c(0.50)) & .x <= quantile(.x, c(0.75)) ~ "med-high",
                                .x > quantile(.x, c(0.75)) ~ "high")))
    cont_col <- "none"
  }
  
  ##make traits long format for joining
  t <- t %>%
    #prep categorical traits to be p/a
    pivot_longer(cols = -c(any_of(tax_col_name), all_of(bi_col), any_of(cont_col)), names_to = "trait", values_to = "value") %>%
    mutate(p = 1) %>%
    pivot_wider(names_from = c(trait, value), values_from = p, values_fill = list(p = 0)) %>%
    pivot_longer(cols = -any_of(tax_col_name), names_to = "trait", values_to = "p") %>%
    filter(trait %in% cont_col | p == 1) #filter out the absences if not a contin. variable
  
  #join site x traits
  if(!"taxa" %in% names(t)) { #rename taxa name column for joining
    names(t)[names(t) %in% tax_col_name] <- "taxa"
  }
  
  if(match_spp_across_dat) {
    t <- t[t$taxa %in% o$taxa, ]
    o <- o[o$taxa %in% t$taxa, ]
  } else {
    #check for matching taxa names if not filtering
    if(!all(o$taxa %in% t$taxa)) {
      warning("Trait taxa and site taxa do not match! You could be losing data.")
    }
  }
  
  t.o <- left_join(o, t, relationship = "many-to-many") %>%
    select(-any_of(tax_col_name))
  
  #split df for making wide df in next step
  cont.t.o <- t.o[t.o$trait %in% cont_col,]
  bi.t.o <- t.o[!t.o$trait %in% cont_col,]
  
  if(nrow(cont.t.o) != 0) {
    bi.t.o <- bi.t.o %>% pivot_wider(names_from = trait, values_from = p, values_fn = ~ sum(.x, na.rm = FALSE))
    cont.t.o <- cont.t.o %>% pivot_wider(names_from = trait, values_from = p, values_fn = ~ mean(.x, na.rm = FALSE))
    
    t.o <- left_join(cont.t.o, bi.t.o)
    
  } else { #if no continuous variables, all variables are treated the same
    t.o <- t.o %>% pivot_wider(names_from = trait, values_from = p, values_fn = ~ sum(.x, na.rm = FALSE))
  }
  
  if(nrow(site_dat) != nrow(t.o)) { message("Oops! Something went wrong. Original # of sites and # of sites in new matrix are not the same.") }
  if(!any(c("presence/absence", "abundance") %in% matrix_type)) { stop("Unsupported matrix type specified. Check spelling.") }
  
  out_list <- vector(mode = "list", length = 2)
  
  if("presence/absence" %in% matrix_type) { #get presence/absence matrix (any species with modality = present)
    pa <- t.o %>%
      select(-tot_tax) %>%
      mutate(across(-c(any_of(ignore_cols), any_of(cont_col)),
                    ~ ifelse(is.na(.x), 0, 1)))
    out_list[[1]] <- pa 
    names(out_list)[1] <- "presence_absence"
  }
  
  if("abundance" %in% matrix_type) {
    abund <- t.o %>%
      mutate(across(.cols = -c(tot_tax, any_of(ignore_cols), any_of(cont_col)), 
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


