##########################################
### Work qPCR data in an automated way ###
##########################################

### WD
setwd("/home/labs/eclab/ijarne/Downloads")

### libraries
library(readr)
library(dplyr)
library(ggplot2)

### Constants
raw_data_file <- "./final_excel_for_L.xlsx"

### pick the data
raw_data <- openxlsx::read.xlsx(xlsxFile = raw_data_file)

### meta data
raw_data <- raw_data %>%
  dplyr::mutate(cell_type = gsub(pattern = ".[0-9]{1}$", replacement = "", x = sample_name),
                cell_type = gsub(pattern = "\\..*$", replacement = "", x = cell_type),
                sample_group = ifelse(grepl(x = sample_name, pattern = "DN"),"Control","Res")) %>%
  dplyr::relocate(cell_type, sample_group, .after = sample_name) 

### clean raw_data and do SD and mean values
final_data <- raw_data %>% 
  dplyr::filter(!is.na(mean_value)) %>% 
  dplyr::

data_to_plot <- raw_data %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(id = paste0(paste0(cell_type,"_",sample_group,"_",Gene))) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(id, .before = 1) %>%
  dplyr::mutate(cell_type = factor(x = cell_type, levels = c("MM074","MM099","MM0383")),
                sample_group = factor(x = sample_group, levels = c("Control","Res")),
                Gene = factor(x = Gene, levels = unique(raw_data$Gene)))
  
  
# data_to_plot_main <- data_to_plot %>% 
#   dplyr::select(-c(Pos, sample_name, type_of_gene, Cp, good_sample, good_sample, n_good, elevated_value)) %>% 
#   dplyr::filter(!is.na(mean_value)) %>% 
#   dplyr::
# 
# 
# plot_of_genes <- ggplot(data = data_to_plot_main, mapping = aes(x = id, y = mean_value,
#                                                                 colour = sample_group, fill = cell_type)) +
#   geom_col()
# 
# plot_of_genes
# 



