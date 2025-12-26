##########################################
### Work qPCR data in an automated way ###
##########################################

### WD
# setwd("/home/labs/eclab/ijarne/Desktop/") ## work
setwd("/home/ignasi/Documents/7qpcr") ## home

### libraries
library(readr)
library(openxlsx)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

## some useful functions
# expand cols
expand_cols <- function(x) {
  unlist(lapply(strsplit(x, ",")[[1]], function(y) {
    if (grepl("-", y)) {
      r <- as.integer(strsplit(y, "-")[[1]])
      seq(r[1], r[2])
    } else {
      as.integer(y)
    }
  }))
}
# expand cols
expand_rows <- function(x) {
  strsplit(x, ",")[[1]]
}
# expand the dataframes
expand_df <- function(df, sup_col = "") {
  expanded_df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      expanded = list(
        expand_grid(
          row = expand_rows(rows),
          col = expand_cols(cols)
        )
      )
    ) %>%
    tidyr::unnest(expanded) %>%
    dplyr::ungroup() %>%
    dplyr::select(file_name,all_of(sup_col), row, col) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Pos = (paste0(row,col))) %>% 
    dplyr:: ungroup() %>% 
    dplyr::select(file_name,all_of(sup_col),Pos) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(id_of_col = paste0(file_name,"_",Pos)) %>%
    dplyr::ungroup() %>% 
    dplyr::relocate(id_of_col, .before = 1)
  
  return(expanded_df)
  }

### Constants
## data placement
raw_data_file <- c("./raw_data_to_play/20251212_SPOCK1_DN-TRES.txt")
sample_location <- c("./raw_data_to_play/sample_location.xlsx")
gene_location <- c("./raw_data_to_play/gene_location.xlsx")

## control group
control_group <- "DN"

## Housekeeping genes
genes_housekeeping <- c("RHOA")

## value for bad sample and bad value
value_for_bad_replicate <- 35
value_for_bad_value <- 0.5
## permited dispersion at the sample level
permited_dispersion <- 3

### Data importation
## raw txt data
raw_data <- lapply(X = raw_data_file, FUN = function(x){
  element_now <- readr::read_tsv(file = x, skip = 1, locale = readr::locale(decimal_mark = ","))
  element_now <- element_now %>% 
    dplyr::mutate(file_name = gsub(pattern = "\\.txt$", 
                                   replacement = "",
                                   x = basename(x))) %>% 
    dplyr::relocate(file_name, .before = 1)
  return(element_now)
  })
raw_data <- do.call(rbind, raw_data)

## sample locations
sample_location <- lapply(X = sample_location, FUN = function(x) {
  element_now <- openxlsx::read.xlsx(xlsxFile = x)
  return(element_now)
  })
sample_location <- do.call(rbind, sample_location)
                
## gene locations
gene_location <- lapply(X = gene_location, FUN = function(x) {
  element_now <- openxlsx::read.xlsx(xlsxFile = x)
  return(element_now)
})
gene_location <- do.call(rbind, gene_location)

### Genes work
## Genes location
gene_location <- expand_df(df = gene_location, sup_col = "gene")
sample_location <- expand_df(df = sample_location, sup_col = "sample")

#### mutate raw_data
raw_data <- raw_data %>%
  # get the useful data
  dplyr::select(c(file_name,Pos,Cp)) %>%
  # generate the id_of_col column
  dplyr::rowwise() %>% 
  dplyr::mutate(id_of_col = paste0(file_name,"_",Pos)) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(id_of_col, .before = 1) %>% 
  # generate col and row data
  dplyr::rowwise() %>% 
  dplyr::mutate(row_name = strsplit(x = Pos, split = "")[[1]][1],
                col_name = strsplit(x = Pos, split = "")[[1]][2]) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(row_name, col_name, .before = Pos) %>% 
  # Add the gene info
  dplyr::mutate(Gene = gene_location$gene[match(id_of_col, gene_location$id_of_col)]) %>%
  dplyr::mutate(type_of_gene = ifelse(Gene %in% genes_housekeeping,"HK","Normal")) %>%
  # # Add replicate_name
  dplyr::mutate(replicate_name = sample_location$sample[match(id_of_col, sample_location$id_of_col)])

raw_data <- raw_data %>% 
  dplyr::select(-c(id_of_col,row_name, col_name)) %>% 
  dplyr::relocate(replicate_name, .after = Pos) %>% 
  dplyr::relocate(Gene, type_of_gene, .after = replicate_name)

########################################################
####                INICIAR LES DADES               ####
########################################################

##### 1 remove replicate_name NA and add sample_name information
raw_data <- raw_data %>% 
  dplyr::filter(!is.na(replicate_name)) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(sample_name = gsub(pattern = "\\.[^.]*$", replacement = "", x = replicate_name)) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(sample_name, .after = replicate_name)

##### 2 mark replicates with HK35 and drop them
raw_data <- raw_data %>% 
  dplyr::group_by(replicate_name) %>%
  dplyr::mutate(
    good_replicate = ifelse(
      any(type_of_gene == "HK" & Cp >= value_for_bad_replicate, na.rm = TRUE),
      "bad","good")) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(good_replicate == "good")

##### 3 NAs to 35
raw_data <- raw_data %>%
  dplyr::mutate(Cp = ifelse(is.na(Cp),value_for_bad_replicate,Cp))

##### 4 check for bad values
raw_data <- raw_data %>%
  # compare each value to check if they are good or bad compared to each of 
  # the replicate
  dplyr::group_by(replicate_name, Gene) %>%
  dplyr::mutate(
    comparison_results = map_chr(Cp, ~ {
      comps <- abs(.x - Cp) <= value_for_bad_value
      comps <- comps[!is.na(comps)]
      comps <- comps[-which(Cp == .x)[1]]  # remove self-comparison

      return(paste0(ifelse(comps, "good", "bad"), collapse = ","))
    })
  ) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(replicate_name, Gene) %>%
  # check the number of times a replicate has values
  # count number of good and bad samples 
  # determine the number of good to keep a sample
  dplyr::mutate(
    valid_counter = str_count(comparison_results, "good"),
    n_of_samples = n(),
    n_of_min_valid = n_of_samples - 2  # threshold
  ) %>%
  # determine the type of sample group we are in
  dplyr::mutate(
    group_type = case_when(
      all(valid_counter == n_of_samples-1) ~ "all_good",
      all(valid_counter == 0)   ~ "all_bad",
      TRUE ~ "mixed")) %>% 
  # determine if we keep it or not
  mutate(
    keep_this_value = case_when(
      group_type %in% c("all_good", "all_bad") ~ "keep",
      group_type == "mixed" & valid_counter >= n_of_min_valid ~ "keep",
      TRUE ~ "NO keep"
    )) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-c(valid_counter,n_of_samples,
                   n_of_min_valid,group_type))

##### 5 mean of the CT
raw_data <- raw_data %>%
  dplyr::group_by(replicate_name, Gene) %>%
  dplyr::mutate(
    mean_value = mean(Cp[keep_this_value == "keep"], na.rm = TRUE),
    sd_value = sd(Cp[keep_this_value == "keep"], na.rm = TRUE)) %>%
  dplyr::ungroup()

##### 6 deltaCT calculation
raw_data <- raw_data %>%
  dplyr::group_by(replicate_name) %>%
  dplyr::mutate(
    hk_mean_value = mean(mean_value[type_of_gene == "HK"],
                         na.rm = TRUE),
    deltaCT = mean_value-hk_mean_value) %>%
  dplyr::ungroup()

openxlsx::write.xlsx(file = "./results/raw_data_replicate_level.xlsx", x = raw_data)

#### 7 delta deltaCT
sample_data <- raw_data %>% 
  # remove non wanted columns
  dplyr::select(-c(file_name,Pos,replicate_name,Cp,good_replicate,comparison_results,keep_this_value,mean_value,sd_value,hk_mean_value)) %>% 
  # remove HK genes
  dplyr::filter(type_of_gene != "HK") %>% 
  dplyr::select(-type_of_gene) %>% 
  # keep only sample level data so remove duplicates
  dplyr::distinct() %>% 
  # add sample group and cell line info
  dplyr::rowwise() %>% 
  dplyr::mutate(sample_group = ifelse(grepl(pattern = control_group, x = sample_name),"Control","Quest"),
                cell_line = strsplit(x = sample_name, split = "\\.")[[1]][1]) %>%
  dplyr::ungroup() %>% 
  dplyr::relocate(cell_line, sample_group, .after = sample_name) %>% 
  # dp the delta delta CT calculation
  dplyr::group_by(sample_name, Gene) %>% 
  dplyr::mutate(deltaCT_mean = mean(deltaCT)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(cell_line,Gene) %>% 
  dplyr::mutate(delta_delta_CT = deltaCT - mean(deltaCT_mean[sample_group == "Control"])) %>% 
  dplyr::ungroup() %>% 
  # calculate elevated value
  dplyr::mutate(elevated_value = 2^-delta_delta_CT) %>% 
  dplyr::select(-c(deltaCT_mean)) %>% 
  # check outliers of elevated value
  dplyr::group_by(cell_line, sample_group, Gene) %>%
  dplyr::mutate(
    mean_others = map_dbl(seq_along(elevated_value), function(i) {
      others <- elevated_value[-i]
      if (length(others) >= 1) mean(others, na.rm = TRUE) else NA_real_
    }),
    sd_others = map_dbl(seq_along(elevated_value), function(i) {
      others <- elevated_value[-i]
      if (length(others) >= 2) sd(others, na.rm = TRUE) else NA_real_
    }),
    qc_flag = ifelse(
      !is.na(sd_others) &
        abs(elevated_value - mean_others) > permited_dispersion * sd_others,
      "BAD",
      "GOOD"
    )
  ) %>%
  dplyr::ungroup() %>% 
  # calculate the mean and sd values of good delta scores
  dplyr::group_by(cell_line,sample_group,Gene) %>%
  dplyr::mutate(mean_value = mean(elevated_value[qc_flag == "GOOD"]),
                sd_value = sd(elevated_value[qc_flag == "GOOD"])) %>% 
  dplyr::ungroup()

##### 8 Clean for Excel and add deltaCT mean and deltaCT sd
sample_data <- sample_data %>%
  dplyr::group_by(cell_line,sample_group,Gene) %>% 
  mutate(
    row_in_group = row_number(),
    mean_value   = ifelse(row_in_group == 1, mean_value, NA_real_),
    sd_value     = ifelse(row_in_group == 1, sd_value, NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::select(-row_in_group) %>% 
  dplyr::arrange(sample_name,sample_group,Gene)

openxlsx::write.xlsx(file = "./results/raw_sample_data.xlsx", x = sample_data)

##### 9 do the zscore calculations of deltaCT
clean_data_to_plot_hmaps <- sample_data %>% 
  # remove bad data and not useful columns
  dplyr::select(-qc_flag) %>% 
  dplyr::filter(!is.na(mean_value)) %>% 
  dplyr::select(sample_name, cell_line,Gene, sample_group, mean_value) %>%
  # calculate the zscore
  dplyr::group_by(cell_line,Gene) %>% 
  dplyr::mutate(z_scored_deltaCT = {
    mean_val <- mean(mean_value)
    sd_val <- sd(mean_value)
    (mean_value-mean_val/sd_val)}) %>% 
  dplyr::ungroup()

openxlsx::write.xlsx(file = "./results/to_heatmaps.xlsx", x = clean_data_to_plot_hmaps)

##### 10 for the barplots
clean_data_to_plot_barplots <- sample_data %>% 
  # remove bad data
  dplyr::filter(qc_flag == "GOOD") %>%
  dplyr::select(-qc_flag) %>% 
  dplyr::ungroup()
  
