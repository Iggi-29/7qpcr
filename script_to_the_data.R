##########################################
### Work qPCR data in an automated way ###
##########################################

### WD
# setwd("/home/labs/eclab/ijarne/Desktop/") ## work
setwd("/home/ignasi/Documents/7qpcr") ## home

### Check and if the  results folder doesn't exist, generate it
if (!dir.exists("./results/")) {
  dir.create(path = "./results/")
}

### libraries
library(readr)
library(openxlsx)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

### functions importation
source("./functions/functions.R")

### Constants
## data placement
raw_data_file <- c("./raw_data/20251212_SPOCK1_DN-TRES.txt")
sample_location <- c("./raw_data/sample_location.xlsx")
gene_location <- c("./raw_data/gene_location.xlsx")
meta_data_location <- "./raw_data/meta_data.xlsx"

## control group
control_group <- "DN"

## Housekeeping genes
genes_housekeeping <- c("RHOA")

## value for bad replicate and bad value
value_for_bad_replicate <- 35
value_for_bad_value <- 0.5

## perform the tests
perform_test <- TRUE

## permited dispersion at the sample level
permited_dispersion <- 2.5

## for the test or wilcoxon test 
test_type <- "ttest"

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
## meta_data
meta_data <- openxlsx::read.xlsx(xlsxFile = meta_data_location)

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

colnames1 <- colnames(raw_data)

########################################################
####                INICIAR LES DADES               ####
########################################################

##### 1 remove replicate_name NA and add meta_data information
raw_data <- raw_data %>% 
  dplyr::filter(!is.na(replicate_name)) 

raw_data <- merge(x = raw_data, y = meta_data,
                  by = "replicate_name")
raw_data <- raw_data %>% 
  dplyr::relocate(colnames(meta_data[2:ncol(meta_data)]), .after = replicate_name) %>% 
  dplyr::relocate(file_name, .before = 1)  

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

openxlsx::write.xlsx(file = "./results/raw_data_sample_location.xlsx",
                     x = raw_data)

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

##### 5 elevated value calculation
raw_data <- raw_data %>% 
  dplyr::mutate(elevated_value = 0.5^Cp)

##### 6 elevated value mean
raw_data <- raw_data %>% 
  dplyr:: group_by(replicate_name,Gene) %>% 
  dplyr::mutate(elevated_value_mean = mean(elevated_value[keep_this_value == "keep"])) %>% 
  dplyr::ungroup()

##### 7 normalized expression
raw_data <- raw_data %>% 
  dplyr::group_by(replicate_name) %>% 
  dplyr::mutate(normalized_expression = {
    mean_hk <- mean(elevated_value_mean[keep_this_value == "keep" & type_of_gene == "HK"])
    elevated_value_mean / mean_hk
  }) %>% 
  dplyr::ungroup() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(log2_expression_value = log2(normalized_expression)) %>% 
  dplyr::ungroup()

##### 8 remove repetitive info at the replicate level
raw_data <- raw_data %>% 
  dplyr::group_by(replicate_name, Gene) %>% 
  dplyr::mutate(elevated_value_mean = ifelse(row_number() == 1,elevated_value_mean,NA),
                normalized_expression = ifelse(row_number() == 1,normalized_expression,NA),
                log2_expression_value = ifelse(row_number() == 1,log2_expression_value,NA)) %>%
  dplyr::ungroup()

##### 9 detect outliers at the sample-level
bio_rep_data <- raw_data %>% 
  dplyr::select(c(
    file_name,Pos,sample_name,zscoring_group,testing_group,
    Gene,type_of_gene,
    elevated_value_mean,normalized_expression,log2_expression_value)) %>% 
  dplyr::filter(!is.na(elevated_value_mean))

bio_rep_data <- bio_rep_data %>%
  dplyr::group_by(sample_name, Gene) %>%
  dplyr::mutate(
    mean_others = map_dbl(seq_along(normalized_expression), function(i) {
      others <- normalized_expression[-i]
      if (length(others) >= 1) mean(others, na.rm = TRUE) else NA_real_
    }),
    sd_others = map_dbl(seq_along(normalized_expression), function(i) {
      others <- normalized_expression[-i]
      if (length(others) >= 2) sd(others, na.rm = TRUE) else NA_real_
    }),
    qc_flag = ifelse(
      !is.na(sd_others) &
        abs(normalized_expression - mean_others) > permited_dispersion * sd_others,
      "BAD",
      "GOOD"
    )
  ) %>%
  dplyr::ungroup()

bio_rep_data_red <- bio_rep_data %>% 
  dplyr::group_by(sample_name, Gene) %>% 
  dplyr:: mutate(normalized_expression_mean = mean(normalized_expression),
                 log2_expression_value_mean = mean(log2_expression_value)) %>% 
  dplyr::mutate(normalized_expression_sd = sd(normalized_expression),
                log2_expression_value_sd = sd(log2_expression_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(file_name,Pos,sample_name,zscoring_group,testing_group,
                  Gene, type_of_gene,
                  normalized_expression_mean, log2_expression_value_mean,
                  normalized_expression_sd, log2_expression_value_sd)) %>% 
  dplyr::distinct()

##### 10 reconstruct the dataset
ttest_ready_data <- bio_rep_data

bio_rep_data <- bio_rep_data %>% 
  dplyr::select(c(file_name,zscoring_group,testing_group,Pos,Gene,mean_others,sd_others,qc_flag))

bio_rep_data_red <- bio_rep_data_red %>% 
  dplyr::group_by(sample_name, Gene) %>% 
  dplyr::mutate(normalized_expression_mean = ifelse(row_number() == 1,normalized_expression_mean,NA),
                normalized_expression_sd = ifelse(row_number() == 1, normalized_expression_sd,NA),
                log2_expression_value_mean = ifelse(row_number() == 1,log2_expression_value_mean,NA),
                log2_expression_value_sd = ifelse(row_number() == 1,log2_expression_value_sd,NA)) %>%
  dplyr::ungroup() %>% 
  dplyr::select(c(file_name,Pos,zscoring_group,testing_group,Gene,
                  normalized_expression_mean,normalized_expression_sd,
                  log2_expression_value_mean,log2_expression_value_sd)) %>% 
  dplyr::filter(!is.na(normalized_expression_mean)) %>% 
  dplyr::group_by(zscoring_group,Gene) %>% 
  dplyr::mutate(zscored_value = {
    mean_val = mean(log2_expression_value_mean, na.rm = TRUE)
    sd_val = sd(log2_expression_value_mean, na.rm = TRUE)
    (log2_expression_value_mean - mean_val)/sd_val
  }) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(file_name,Pos,normalized_expression_mean,normalized_expression_sd,
                  log2_expression_value_mean,log2_expression_value_sd,zscored_value))

bio_rep_data <- bio_rep_data %>% 
  dplyr::select(c(file_name,Pos,mean_others,sd_others,qc_flag))

raw_data <- merge(
  x = raw_data, y = bio_rep_data,
  by = c("file_name","Pos"), all.x = TRUE)
raw_data <- merge(
  x = raw_data, y = bio_rep_data_red,
  by = c("file_name","Pos"), all.x = TRUE)

raw_data <- raw_data %>% 
  dplyr::arrange(replicate_name, Gene, is.na(elevated_value_mean))

##### 11 export the final data
### Generate the wb structure
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "Sheet1")
openxlsx::writeData(wb = wb, sheet = "Sheet1", x = raw_data)
### Headers
## annotation cols
header_style <- createStyle(
  fgFill = "white",   # background color
  textDecoration = "bold",
  halign = "left",
  border = c("Bottom","Right","Left","Top"))
addStyle(
  wb, sheet = "Sheet1",
  style = header_style,
  rows = 1, 
  cols = grep(x = colnames(raw_data),
              pattern = paste0(x = c(colnames1, colnames(meta_data)),
                               collapse = "|")),
  gridExpand = TRUE)

### Headers
## annotation cols
header_style <- createStyle(
  fgFill = "white",   # background color
  textDecoration = "bold",
  halign = "left",
  border = c("Bottom","Right","Left","Top"))
addStyle(
  wb, sheet = "Sheet1",
  style = header_style,
  rows = 1, 
  cols = grep(x = colnames(raw_data),
              pattern = paste0(x = c(colnames1, colnames(meta_data)),
                               collapse = "|")),
  gridExpand = TRUE)

## replicate cols
header_style <- createStyle(
  fgFill = "lightgreen",
  textDecoration = "bold",
  halign = "left",
  border = c("Bottom","Right","Left","Top"))
addStyle(
  wb, sheet = "Sheet1",
  style = header_style,
  rows = 1, 
  cols = grep(x = colnames(raw_data),
              pattern = paste0(x = c("good_replicate","comparison_results",
                                     "keep_this_value","elevated_value"),
                               collapse = "|")),
  gridExpand = TRUE)
## replicate goruped cols
header_style <- createStyle(
  fgFill = "lightgreen",
  textDecoration = "bold",
  halign = "left",
  border = c("Bottom","Right","Left","Top"))
addStyle(
  wb, sheet = "Sheet1",
  style = header_style,
  rows = 1, 
  cols = grep(x = colnames(raw_data),
              pattern = paste0(x = c("normalized_expression","log2_expression_value",
                                     "mean_others","sd_others","qc_flag"),
                               collapse = "|")),
  gridExpand = TRUE)
## sample goruped cols
header_style <- createStyle(
  fgFill = "lightblue",
  textDecoration = "bold",
  halign = "left",
  border = c("Bottom","Right","Left","Top"))
addStyle(
  wb, sheet = "Sheet1",
  style = header_style,
  rows = 1, 
  cols = grep(x = colnames(raw_data),
              pattern = paste0(x = c("normalized_expression_mean","normalized_expression_sd",
                                     "log2_expression_value_mean","log2_expression_value_sd",
                                     "zscored_value"),
                               collapse = "|")),
  gridExpand = TRUE)


saveWorkbook(wb, "./results/final_data.xlsx", overwrite = TRUE)

#### 12 do the t-est or wald test
if (perform_test)  {
  ttest_ready_data <- ttest_ready_data %>%
    dplyr::select(sample_name,zscoring_group,testing_group,Gene,type_of_gene,normalized_expression) %>%
    dplyr::filter(type_of_gene != "HK") %>%
    dplyr::group_by(zscoring_group,Gene) %>%
    dplyr::summarise(pvalue = run_test(
      df = pick(everything()),
      value_col = "normalized_expression",
      group_col = "testing_group",
      test_type = test_type)) %>%
    dplyr::ungroup()
  
  openxlsx::write.xlsx(file = "./results/statistical_results.xlsx", 
                       x = ttest_ready_data)
  }

