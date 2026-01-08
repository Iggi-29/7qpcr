 
### Read the data
heatmap_data <- openxlsx::read.xlsx(xlsxFile = heatmap_data_file)
gene_annotation <- openxlsx::read.xlsx(xlsxFile = gene_annotation_file)

## some useful functions
# expand cols

### Generate general objects for heatmap
# color palette for the heatmap
colours_of_allota <- c(
  "#FFF5F5","#F7E3E3",
  "#EEC6C6","#E6AAAA",
  "#DD8E8E","#D57171",
  "#CC5555","#C43939",
  "#BB1C1C","#B30000")

# color of genes

# order of the data

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


