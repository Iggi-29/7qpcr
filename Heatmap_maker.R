################################################
## Zscore qPCR results and generate a heatmap ##
################################################

### libraries and WD
library(dplyr)
library(tidyr)
library(readxl)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# setwd("/Users/ltorres/Desktop/")
setwd("/home/ignasi/Desktop/heatmaps/")

# color palette for the heatmap
colours_of_allota <- c("#FFF5F5","#F7DADA", 
                       "#EEBFBF","#E6A3A3",
                       "#DD8888","#D56D6D", 
                       "#CC5252","#C43636",
                       "#BB1B1B","#B30000")


colours_of_allota <- c(
  "#FFF5F5",
  "#F7E3E3",
  "#EEC6C6",
  "#E6AAAA",
  "#DD8E8E",
  "#D57171",
  "#CC5555",
  "#C43939",
  "#BB1C1C",
  "#B30000"
)
### read the data
sheet_names <- readxl::excel_sheets(path = "./data/Z-Score DTPs.xlsx")
raw_data <- lapply(X = sheet_names, FUN = function(x) {
  readxl::read_xlsx(path = "./data/Z-Score DTPs.xlsx", sheet = x)})
names(raw_data) <- sheet_names

# clean colnames just in case ... this shall be removed maybe .. but it is not
# necessary
for (i in 1:length(raw_data)) {
  names_now <- names(raw_data)[i]
  raw_data_now <- raw_data[[i]]

  colnames(raw_data_now) <- gsub(pattern = "\\.\\.\\.",
                                 replacement = "_R_",
                                 x = colnames(raw_data_now))
  names(raw_data)[i] <- names_now
  raw_data[[i]] <- raw_data_now}


zscored_list <- list()
## zscore in a very easy way
for (i in 1:length(raw_data)) {
  names_now <- names(raw_data)[i]
  raw_data_now <- raw_data[[i]]
  
  # to long format
  raw_data_now_long <- tidyr::pivot_longer(
    data = raw_data_now,
    cols = colnames(raw_data_now)[2:ncol(raw_data_now)],
    values_to = "qvalue", names_to = "sample_name")
  
  # remove NAs
  raw_data_now_long <- raw_data_now_long %>% 
    dplyr::filter(!is.na(qvalue))
  
  # detect min value
  raw_data_now_long <- raw_data_now_long %>% 
    dplyr::mutate(qvalue = ifelse(qvalue == min(raw_data_now_long$qvalue),NA,qvalue))
  
  # replcates into samples
  raw_data_now_long <- raw_data_now_long %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(sample_name = gsub(pattern = "_.*", replacement = "", x = sample_name)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(GENES,sample_name) %>% 
    dplyr::mutate(mean_qvalue = mean(qvalue)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-qvalue) %>% 
    dplyr::rename(qvalue = mean_qvalue) %>% 
    dplyr::distinct()
  
  # zscore values
  zscored_values <- raw_data_now_long %>%
    # calculate the mean and sd values
    dplyr::group_by(GENES) %>%
    dplyr::mutate(mean_val = mean(qvalue, na.rm = T),
                  sd_val = sd(qvalue, na.rm = T)) %>%
    dplyr::ungroup() %>%
    # calculate the zscore
    dplyr::rowwise() %>% 
    dplyr::mutate(zscore_value = (qvalue-mean_val)/sd_val) %>% 
    dplyr::ungroup()
  
  zscored_list[[i]] <- zscored_values
  names(zscored_list)[i] <- names_now
  }

## do the heatmaps
for (i in 1:length(zscored_list)) {
  # data now
  zscored_list_now <- zscored_list[[i]]
  names_now <- names(zscored_list)[i]
  
  # pick the data
  normal_data <- zscored_list_now %>% 
    dplyr::select(GENES,sample_name,qvalue) %>% 
    tidyr::pivot_wider(
      names_from = sample_name,
      values_from = qvalue)
  normal_data <- as.data.frame(normal_data)
  rownames(normal_data) <- normal_data$GENES
  normal_data <- normal_data[,-(grep(pattern = "GENES", x = colnames(normal_data))), drop = FALSE]
  
  zscored_data <- zscored_list_now %>% 
    dplyr::select(GENES,sample_name,zscore_value) %>% 
    tidyr::pivot_wider(
      names_from = sample_name,
      values_from = zscore_value) 
  zscored_data <- as.data.frame(zscored_data)
  rownames(zscored_data) <- zscored_data$GENES
  zscored_data <- zscored_data[,-(grep(pattern = "GENES", x = colnames(zscored_data))), drop = FALSE]
  
  # heatmaps
  normal_hm <- ComplexHeatmap::Heatmap(
    matrix = as.matrix(normal_data),
    cluster_columns = FALSE, cluster_rows = FALSE, 
    na_col = "white", column_names_rot = 90,
    col = colours_of_allota, 
    row_names_side = "left",
    name = "Raw\nexpression value", 
    column_title = names_now)

  zscored_hm <- ComplexHeatmap::Heatmap(
    matrix = as.matrix(zscored_data),
    cluster_columns = FALSE, cluster_rows = FALSE, 
    na_col = "white", column_names_rot = 90,
    col = colours_of_allota,
    row_names_side = "left",
    name = "Z-scored\nexpression value", 
    column_title = names_now)
  
  # save the heatmaps - PDF
  pdf(file = paste0("./plots/",names_now,"_normal_data.pdf"), width = 6, height = 8)
  normal_hm <- ComplexHeatmap::draw(normal_hm)
  dev.off()

  pdf(file = paste0("./plots/",names_now,"_zscored_data.pdf"), width = 6, height = 8)
  zscored_hm <- ComplexHeatmap::draw(zscored_hm)
  dev.off()
  
  # save the heatmaps - png
  png(filename = paste0("./plots/",names_now,"_normal_data.png"), width = 6, height = 8, units = "in", res = 300)
  normal_hm <- ComplexHeatmap::draw(normal_hm)
  dev.off()
  
  png(filename = paste0("./plots/",names_now,"_zscored_data.png"), width = 6, height = 8, units = "in", res = 300)
  zscored_hm <- ComplexHeatmap::draw(zscored_hm)
  dev.off()
  
  }


