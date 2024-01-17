### COMPARE GF MODEL OUTPUTS ###
#make summary figures for comparing variable types, taxa, gage types etc.

library(tidyverse); library(gradientForest)
options(readr.show_col_types = FALSE)

PATH <- getwd()

gf.list <- list.files(paste0(PATH, "/10_GFOutput/2023_11_14"), ".rds", full.names = TRUE)

#set colors for consistency ----
source(paste0(PATH, "/Scripts/XX_colors.R"))

gf_file <- gf.list[[1]] #testing

# variable importance (weighted R2) ----
v.imp <- lapply(gf.list, function(gf_file) { #get raw weighted R2 data
  
  gf.out <- readRDS(gf_file)
  imp.vars <- data.frame(importance(gf.out)) %>%
    rownames_to_column("env_var")
  names(imp.vars)[[2]] <- "weighted_r2"
  
  gf_name <- sub(".rds", "", sub(".*/", "", gf_file))
  
  imp.vars <- imp.vars %>%
    mutate(taxa = case_when(grepl("bug", gf_name) ~ "bug",
                            grepl("fish", gf_name) ~ "fish"),
           flow = case_when(grepl("GW", gf_name) ~ "GW",
                            grepl("Int", gf_name) ~ "Int",
                            grepl("RO", gf_name) ~ "RO"),
           gage = case_when(grepl("non-refgage", gf_name) ~ "non-ref",
                            grepl("allgage", gf_name) ~ "all",
                            T ~ "ref"),
           var_grp = case_when(grepl("HIT", gf_name) ~ "HIT",
                               grepl("HyALT", gf_name) ~ "HyALT",
                               grepl("LULC", gf_name) ~ "LULC",
                               T ~ "all"))
  
  return(imp.vars)
  
})

v.imp <- bind_rows(v.imp)

write_csv(v.imp, paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/variable_importance_R2_combined.csv"))

##plot VI (weighted R2) ----
v.imp <- read_csv(paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/variable_importance_R2_combined.csv"))

#get variable types for colors
v.imp <- v.imp %>% 
  left_join(., 
            v.imp %>%
              select(env_var, var_grp) %>%
              filter(var_grp != "all") %>%
              rename(var_type_lrg = var_grp) %>%
              distinct()
            ) %>%
  mutate(var_type_sml = case_when(var_type_lrg == "HIT" ~ "hydrology",
                                  grepl("BFIWs", env_var) ~ "hydrology",
                                  grepl("ssn|cv", env_var) ~ "stream temp",
                                  grepl("TmaxWs|TminWs|TmeanWs|PrecipWs", env_var) ~ "climate",
                                  grepl("pn", env_var) ~ "hydrologic alteration",
                                  grepl("Dam|Rd|Pct|NABD|NWs", env_var) ~ "land use",
                                  grepl("Elev|Area|site_", env_var) ~ "spatial",
                                  grepl("Clay|Sand|Om|Perm|WtDep", env_var) ~ "soil",
                                  grepl("Rck|HydrlCond|P2O5|K2O|Kffact", env_var) ~ "geology/lithology",
                                  T ~ "?"))

taxa.vect <- unique(v.imp$taxa)
var.vect <- unique(v.imp$var_grp)

for(sel.taxa in taxa.vect) {
  for(v.group in var.vect) {
    
    #get order by weighted R2 within each flow category
    fl.lvl <- v.imp %>%
      filter(taxa == sel.taxa & var_grp == v.group & gage == "all") %>%
      group_by(flow) %>%
      do(data_frame(al=levels(reorder(interaction(.$flow, .$env_var, drop=TRUE), .$weighted_r2)))) %>%
      pull(al)
    #plot
    plot.df <- v.imp %>%
      filter(taxa == sel.taxa & var_grp == v.group & gage == "all") %>%
      arrange(desc(weighted_r2)) %>%
      group_by(flow) %>%
      slice_max(weighted_r2, n = 25) %>%
      mutate(al = factor(interaction(flow, env_var), levels=fl.lvl))
    p1 <- ggplot(data = plot.df, aes(x = weighted_r2, y = al, fill = var_type_sml)) +
      geom_col() +
      facet_grid(vars(flow), scales = "free_y") +
      scale_x_continuous(expand = c(0.001,0), limits = c(0, max(plot.df$weighted_r2, na.rm = TRUE)+0.0005)) +
      scale_y_discrete(breaks = fl.lvl, labels = sub("^[^.]+\\.", "", fl.lvl)) +
      scale_fill_manual(values = c.pal, name = "Variable Type") +
      theme_bw() + ggtitle(paste0(str_to_title(sel.taxa), ": ", v.group, " variables")) + ylab("")
    
    ggsave(paste0(PATH, "/10_GFOutput/2023_11_14/figures/", sel.taxa, "_variable_importance_", v.group, "vars.png"), plot = p1, width = 8, height = 9)
    
  }
}

rm(p1, plot.df, c.pal, fl.lvl, sel.taxa, taxa.vect, v.group, var.vect)

#cumulative importance ----
ci <- lapply(gf.list, function(gf_file) { #get raw CI data for all models
  
  gf.out <- readRDS(gf_file)
  vars <- names(gf.out$X)
  
  ci.max <- data.frame()
  for(v in vars) {
    
    ci1 <- data.frame(env_var = v, 
                      cum_imp = max(cumimp(gf.out, predictor = v, "Overall")$y))
    
    ci.max <- bind_rows(ci.max, ci1)
    
  }
  
  gf_name <- sub(".rds", "", sub(".*/", "", gf_file))
  
  ci.max <- ci.max %>%
    mutate(taxa = case_when(grepl("bug", gf_name) ~ "bug",
                            grepl("fish", gf_name) ~ "fish"),
           flow = case_when(grepl("GW", gf_name) ~ "GW",
                            grepl("Int", gf_name) ~ "Int",
                            grepl("RO", gf_name) ~ "RO"),
           gage = case_when(grepl("non-refgage", gf_name) ~ "non-ref",
                            grepl("allgage", gf_name) ~ "all",
                            T ~ "ref"),
           var_grp = case_when(grepl("HIT", gf_name) ~ "HIT",
                               grepl("HyALT", gf_name) ~ "HyALT",
                               grepl("LULC", gf_name) ~ "LULC",
                               T ~ "all"))
  
  return(ci.max)
  
})

ci <- bind_rows(ci)

#get variable types for grouping
ci <- ci %>% 
  left_join(., 
            ci %>%
              select(env_var, var_grp) %>%
              filter(var_grp != "all") %>%
              rename(var_type_lrg = var_grp) %>%
              distinct()
  ) %>%
  mutate(var_type_sml = case_when(var_type_lrg == "HIT" ~ "hydrology",
                                  grepl("BFIWs", env_var) ~ "hydrology",
                                  grepl("ssn|cv", env_var) ~ "stream temp",
                                  grepl("TmaxWs|TminWs|TmeanWs|PrecipWs", env_var) ~ "climate",
                                  grepl("pn", env_var) ~ "hydrologic alteration",
                                  grepl("Dam|Rd|Pct|NABD|NWs", env_var) ~ "land use",
                                  grepl("Elev|Area|site_", env_var) ~ "spatial",
                                  grepl("Clay|Sand|Om|Perm|WtDep", env_var) ~ "soil",
                                  grepl("Rck|HydrlCond|P2O5|K2O|Kffact", env_var) ~ "geology/lithology",
                                  T ~ "?"))

write_csv(ci, paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/cumulative_importance_combined.csv"))

##plot CI ----
ci <- read_csv(paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/cumulative_importance_combined.csv"))

##+ all taxa, all var models ----
plot.df <- ci %>%
  filter(var_grp == "all" & gage == "all") %>%
  mutate(var_type_sml = case_when(var_type_sml == "hydrologic alteration" ~ "hydrologic\nalteration", T ~ var_type_sml)) %>%
  group_by(taxa, var_type_sml) %>%
  mutate(mean_value = mean(cum_imp)) %>%
  ungroup() %>%
  arrange(mean_value) %>%
  mutate(flow = factor(flow, levels = c("GW", "RO", "Int")),
         var_type_sml = factor(var_type_sml, levels = unique(var_type_sml), ordered = TRUE)) %>%
  select(-mean_value)

ggplot(data = plot.df, aes(x = var_type_sml, y = cum_imp, fill = flow)) +
  geom_boxplot(position = "dodge", outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.75),
             size = 1,
             color = "black",
             alpha = 0.6) +
  facet_wrap(~ taxa) +
  scale_fill_manual(values = flow.pal, name = "Flow Class") +
  labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
        panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed"))

ggsave(paste0(PATH, "/10_GFOutput/2023_11_14/figures/alltax_cumulative_importance_allvars.png"), width = 16, height = 8)

##+ gage vs. type var models ----
taxa.vect <- unique(ci$taxa)
for(sel.taxa in taxa.vect) {
  
  plot.df <- ci %>%
    filter(taxa == sel.taxa & var_grp == "all" & gage != "all") %>%
    mutate(var_type_sml = case_when(var_type_sml == "hydrologic alteration" ~ "hydrologic\nalteration", T ~ var_type_sml)) %>%
    group_by(var_type_sml) %>%
    mutate(mean_value = mean(cum_imp)) %>%
    ungroup() %>%
    arrange(mean_value) %>%
    mutate(flow = factor(flow, levels = c("GW", "RO", "Int")),
           var_type_sml = factor(var_type_sml, levels = unique(var_type_sml), ordered = TRUE)) %>%
    select(-mean_value)
  
  p2 <- ggplot(data = plot.df, aes(x = var_type_sml, y = cum_imp, fill = gage)) +
    geom_boxplot(position = "dodge", outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.75),
               size = 1,
               color = "black",
               alpha = 0.6) +
    facet_wrap(~ flow) +
    scale_fill_manual(values = gage.pal, name = "Gage Type") +
    labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")")) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
          panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")) +
    ggtitle(paste0(str_to_title(sel.taxa), ": all variables"))
  
  ggsave(paste0(PATH, "/10_GFOutput/2023_11_14/figures/", sel.taxa, "_cumulative_importance_allvars.png"), plot = p2, width = 16, height = 9)
  
}

rm(p2, plot.df, sel.taxa, taxa.vect)

#cumulative importance curves ----

ci.cur <- lapply(gf.list, function(gf_file) { #get raw CI data for all models
  
  gf.out <- readRDS(gf_file)
  vars <- names(gf.out$X)
  
  ci.df <- data.frame()
  for(v in vars) {
    
    xy <- cumimp(gf.out, predictor = v, "Overall")
    
    ci1 <- data.frame(env_var = v, 
                      var_val = xy$x,
                      cum_imp = xy$y)
    
    ci.df <- bind_rows(ci.df, ci1)
    
  }
  
  gf_name <- sub(".rds", "", sub(".*/", "", gf_file))
  
  ci.df <- ci.df %>%
    mutate(taxa = case_when(grepl("bug", gf_name) ~ "bug",
                            grepl("fish", gf_name) ~ "fish"),
           flow = case_when(grepl("GW", gf_name) ~ "GW",
                            grepl("Int", gf_name) ~ "Int",
                            grepl("RO", gf_name) ~ "RO"),
           gage = case_when(grepl("non-refgage", gf_name) ~ "non-ref",
                            grepl("allgage", gf_name) ~ "all",
                            T ~ "ref"),
           var_grp = case_when(grepl("HIT", gf_name) ~ "HIT",
                               grepl("HyALT", gf_name) ~ "HyALT",
                               grepl("LULC", gf_name) ~ "LULC",
                               T ~ "all"))
  
  return(ci.df)
  
})

ci.cur <- bind_rows(ci.cur)

#get variable types for grouping
ci.cur <- ci.cur %>% 
  left_join(., 
            ci.cur %>%
              select(env_var, var_grp) %>%
              filter(var_grp != "all") %>%
              rename(var_type_lrg = var_grp) %>%
              distinct()
  ) %>%
  mutate(var_type_sml = case_when(var_type_lrg == "HIT" ~ "hydrology",
                                  grepl("BFIWs", env_var) ~ "hydrology",
                                  grepl("ssn|cv", env_var) ~ "stream temp",
                                  grepl("TmaxWs|TminWs|TmeanWs|PrecipWs", env_var) ~ "climate",
                                  grepl("pn", env_var) ~ "hydrologic alteration",
                                  grepl("Dam|Rd|Pct|NABD|NWs", env_var) ~ "land use",
                                  grepl("Elev|Area|site_", env_var) ~ "spatial",
                                  grepl("Clay|Sand|Om|Perm|WtDep", env_var) ~ "soil",
                                  grepl("Rck|HydrlCond|P2O5|K2O|Kffact", env_var) ~ "geology/lithology",
                                  T ~ "?"))

write_csv(ci.cur, paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/cumulative_importance_curve_combined.csv"))

##plot CI curve ----
ci.cur <- read_csv(paste0(PATH, "/10_GFOutput/2023_11_14/results_tables/cumulative_importance_curve_combined.csv"))

taxa.vect <- unique(ci.cur$taxa)
var.vect <- unique(ci.cur$var_grp)
gage.vect <- c("all", "all|ref")

# sel.taxa <- "fish"
# v.group <- "HyALT"
# g.group <- "all|ref"

for(sel.taxa in taxa.vect) {
  for(v.group in var.vect) { 
    for(g.group in gage.vect) {
      
      ci.cur.thm <- list(
        labs(y = bquote("Cumulative Importance (" * italic(R[" c"]^2) * ")")),
        theme_bw(),
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", angle = 35, hjust = 0.9, vjust = 1.01),
          panel.grid.major.y = element_line(color = "lightgrey", linetype = "dashed")
        )
      )
      
      plot.df <- ci.cur %>%
        filter(taxa == sel.taxa & var_grp == v.group & grepl(g.group, gage))
      #get top vars
      max.cat <- plot.df %>%
        group_by(flow, env_var) %>%
        slice_max(cum_imp, n = 1) %>%
        ungroup() %>%
        group_by(flow) %>%
        slice_max(cum_imp, n = 4) %>%
        ungroup() %>%
        mutate(imp_for_GW = ifelse(flow == "GW", TRUE, NA),
               imp_for_Int = ifelse(flow == "Int", TRUE, NA),
               imp_for_RO = ifelse(flow == "RO", TRUE, NA)) %>%
        select(env_var, contains("imp_for")) %>%
        group_by(env_var) %>%
        summarise_all(~ any(. == TRUE, na.rm = TRUE)) %>%
        ungroup() %>%
        pivot_longer(cols = contains("imp"), names_to = "flow", values_to = "imp_for") %>%
        mutate(flow = sub("imp_for_", "", flow),
               imp_for = ifelse(imp_for == TRUE, TRUE, NA))
      plot.df <- plot.df %>%
        filter(env_var %in% max.cat$env_var) %>%
        left_join(., max.cat) %>%
        mutate(imp_for = ifelse(is.na(imp_for), "not", "imp. variable"))
      
      if(g.group == "all") {
        
        p1 <- ggplot(data = plot.df, aes(x = var_val, y = cum_imp)) +
          geom_line(aes(color = flow, linetype = imp_for)) +
          scale_color_manual(values = flow.pal, name = "Flow Class") +
          scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
          facet_wrap(~env_var, scales = "free_x") +
          ci.cur.thm +
          ggtitle(paste0(str_to_title(sel.taxa), ": ", v.group, " variables, all gage model"))
        # geom_point(data = max.cat, aes(x = 1, y = 0.020, color = flow))
        
        ggsave(paste0(PATH, "/10_GFOutput/2023_11_14/figures/", sel.taxa, "_cumulative_importance_curve_", v.group, "vars_allgagemods.png"), plot = p1, width = 9, height = 9)
        
      } else {
        
        #ref vs non ref plots
        p1 <- ggplot(data = plot.df, aes(x = var_val, y = cum_imp)) +
          geom_line(aes(color = gage, linetype = imp_for)) +
          scale_color_manual(values = gage.pal, name = "Gage Type") +
          scale_linetype_manual(values = c("imp. variable" = "solid", "not" = "dashed"), name = "") +
          facet_grid(flow~env_var, scales = "free") +
          # geom_text(data = filter(max.cat, imp_for), aes(x = -Inf, y = Inf, label = "*"), hjust = -0.05, vjust = 1.1, size = 18, inherit.aes = FALSE) +
          ci.cur.thm +
          ggtitle(paste0(str_to_title(sel.taxa), ": ", v.group, " variables, gage subset comparisons"))
        
        ggsave(paste0(PATH, "/10_GFOutput/2023_11_14/figures/", sel.taxa, "_cumulative_importance_curve_", v.group, "vars_gagecomp.png"), plot = p1, width = 15, height = 8)
        
        }
      }
    } 
  }

rm(sel.taxa, v.group, g.group, ci.cur.thm, max.cat, p1, plot.df, taxa.vect, var.vect, gage.vect)



#SCRATCH:
# ##note: because of the way the function is set up, need to run the following lines before running split density plot
# gf.out$call$compact <- TRUE
# gf.out$call$nbin <- 201
# plot(gf.out, plot.type = "S", imp.vars = most_important.gf.GW_lump,
#      leg.posn = "topright", cex.legend = 0.9, cex.axis = 0.9,
#      cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
# # 
# # 
# plot(gf.out, plot.type = "Split.Density")
