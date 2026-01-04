##########################################
### Work qPCR data in an automated way ###
##########################################

### WD
# setwd("/home/labs/eclab/ijarne/Desktop/") ## work
setwd("/home/ignasi/Documents/7qpcr") ## home

### Check and if the  results folder doesn't exist, generate it
if (!dir.exists("./plots/barplots")) {
  dir.create(path = "./plots/barplots", recursive = TRUE)
}

### libraries
library(readr)
library(openxlsx)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(ggbreak)

### functions
source("./functions/functions.R")

### Constants
## data placement
raw_data_file <- "./results/final_data.xlsx"
ttest_result_file <- "./results/statistical_results.xlsx"

gene_annotation_file <- "./raw_data/gene_annotation.xlsx"
sample_annotation_file <- "./raw_data/sample_annotation.xlsx" 
colour_of_groups_file <- "./raw_data/colour_of_groups.xlsx"

## control group
control_group <- "DN"

## Housekeeping genes
genes_housekeeping <- c("RHOA")

## jump
jump <- 2

## for the test or wilcoxon test 
test_type <- "ttest"

## Recalculate
recalcultate <- TRUE

### Data importation
## raw_data
raw_data <- openxlsx::read.xlsx(xlsxFile = raw_data_file)

## ttest results
ttest_result <- openxlsx::read.xlsx(xlsxFile = ttest_result_file)

## gene annotation
gene_annotation <- openxlsx::read.xlsx(xlsxFile =  gene_annotation_file)

## sample annotation
sample_annotation <- openxlsx::read.xlsx(xlsxFile = sample_annotation_file)

## colour of groups
colour_of_groups <- openxlsx::read.xlsx(xlsxFile = colour_of_groups_file)
colour_of_groups <- setNames(colour_of_groups$colour, colour_of_groups$testing_group)

### Recalculate if required
if (recalcultate) {
  source("./functions/epp.R")
}
if ("mean_others" %in% colnames(raw_data)) {
  raw_data <- raw_data %>%
    dplyr::select(-c(mean_others,sd_others,comparison_results))
}

### Prepare the data
data_to_plot <- raw_data %>% 
  ## only good columns
  dplyr::select(-c(file_name,Pos,replicate_name,Cp,good_replicate,keep_this_value,elevated_value,zscored_value)) %>% 
  dplyr::select(-c(elevated_value_mean,log2_expression_value,log2_expression_value_mean,log2_expression_value_sd)) %>% 
  ## remove HK genes
  dplyr::filter(type_of_gene !=  "HK") %>% 
  dplyr::select(-c(type_of_gene)) %>% 
  ## remove NA columns
  dplyr::filter(!is.na(normalized_expression))

## do the plots
for (i in 1:length(unique(data_to_plot$zscoring_group)))  {
# for (i in 1:1) {
  ## data now
  group_now <- unique(data_to_plot$zscoring_group)[i]
  data_now <- data_to_plot %>% 
    dplyr::filter(zscoring_group == group_now)
  ttest_result_now <- ttest_result %>% 
    dplyr::filter(zscoring_group == group_now) %>% 
    dplyr::mutate(p_label = paste0("p = ",signif(pvalue,2)))
  sample_annotation_now <- sample_annotation %>% 
    dplyr::filter(zscoring_group == group_now)
  ## order the data
  data_now <- data_now %>% 
    dplyr::arrange(Gene,testing_group) %>%
    dplyr::mutate(ordering = paste0(sample_name,"_",testing_group,"_",Gene)) %>% 
    dplyr::mutate(ordering = factor(x = ordering, levels = unique(ordering))) %>% 
    dplyr::ungroup() %>% 
    dplyr::relocate(ordering, .before = 1)
  ## generate data_bars
  data_bars <- data_now %>% 
    dplyr::filter(!is.na(normalized_expression_mean))
  # y_positions <- data_bars %>%
  #   dplyr::group_by(Gene) %>%
  #   dplyr::summarise(
  #     y_pos = max(normalized_expression_mean + normalized_expression_sd, na.rm = TRUE) * 3,
  #     .groups = "drop")
  y_positions <- max(data_now$normalized_expression_mean, na.rm = TRUE)
  
  ttest_result_now <- ttest_result_now %>%
    dplyr::mutate(y_pos = y_positions+ (y_positions*0.1)) %>% 
    # dplyr::left_join(y_positions, by = "Gene") %>%
    dplyr::mutate(x_pos = median(x = c(1,length(colour_of_groups))))
  
  ### generate the breaks!
  all_points <- data.frame(
    points_of_plot = c(ttest_result_now$y_pos,
                       data_now$normalized_expression)) 
  all_points <- all_points %>% 
    dplyr::arrange(points_of_plot)
  brekk_up <- detect_frontiers(x = all_points$points_of_plot, magnitude_jump = jump)
  brekk_down <- brekk_up -1
  all_points <- all_points[rownames(all_points) %in% c(brekk_down,brekk_up),, drop = FALSE]
  all_points$id <- as.numeric(rownames(all_points))
  all_points <- all_points %>% 
    dplyr::mutate(up_down = ifelse(id %in% brekk_up,"UP",
                                   ifelse(id %in% brekk_down,"DOWN",NA))) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(real_break = ifelse(up_down == "UP",points_of_plot-(points_of_plot/10),
                                      ifelse(up_down == "DOWN",points_of_plot+(points_of_plot/10)))) %>% 
    dplyr::ungroup()
  
  ## do the plot
  bar_plot_cool <- ggplot(data = data_bars, mapping = aes(x = testing_group, y = normalized_expression_mean)) +
    geom_errorbar(mapping = aes(ymin = normalized_expression_mean,
                                ymax = normalized_expression_mean + normalized_expression_sd),
                  linewidth = 0.8) +
    geom_col(mapping = aes(colour = testing_group), 
             fill = sample_annotation_now$colour,
             linewidth = 0.8)+
    scale_colour_manual(values = colour_of_groups) +
    geom_point(data = data_now %>% 
                 dplyr::filter(qc_flag == "GOOD"),
               mapping = aes(x = testing_group, y = normalized_expression),
               shape = 16, size = 3, colour = "black") +
    # geom_point(data = data_now %>%
    #              dplyr::filter(qc_flag == "BAD"),
    #            mapping = aes(x = testing_group, y = normalized_expression),
    #            shape = 16, size = 3, colour = "blue") +
    geom_text(data = ttest_result_now,
              aes(x = x_pos, y = y_pos, label = p_label),
              hjust = 0.5,vjust = 0,size = 4,inherit.aes = FALSE) +
    facet_grid(. ~ Gene, switch = "x") +
    xlab("Gene") +
    ylab("Relative Gene Expression") +
    labs(title = group_now, colour = "Condition") +
    theme_classic() +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  if (nrow(all_points) > 1) {
    bar_plot_cool <- bar_plot_cool +
      scale_y_break(breaks = all_points$real_break, scales = "free")
  }
    ggsave(filename = paste0("./plots/barplots/",group_now,"_all_samples.png"), plot = bar_plot_cool, width = 10, height = 10)
    ggsave(filename = paste0("./plots/barplots/",group_now,"_all_samples.pdf"), plot = bar_plot_cool, width = 10, height = 10)
    
  }  






