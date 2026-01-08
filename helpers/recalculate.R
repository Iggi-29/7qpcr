### Modify the raw_data object to coincide with the new flags
## Select the columns that have information
raw_data2 <- raw_data %>% 
  dplyr::select(file_name,Pos,replicate_name,sample_name,zscoring_group,testing_group,
                Gene,type_of_gene,
                Cp,
                good_replicate,keep_this_value,
                elevated_value,qc_flag)  
  
## remake the information by replicate
raw_data2 <- raw_data2 %>% 
  dplyr::group_by(replicate_name,Gene) %>% 
  dplyr::mutate(elevated_value_mean = mean(elevated_value[keep_this_value == "keep"])) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(replicate_name) %>% 
  dplyr::mutate(normalized_expression = {
    mean_hk = mean(elevated_value_mean[keep_this_value == "keep" & type_of_gene == "HK"])
    elevated_value_mean / mean_hk
  }) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(log2_expression_value = log2(normalized_expression)) %>% 
  dplyr::group_by(replicate_name,Gene) %>% 
  dplyr::mutate(elevated_value_mean = ifelse(row_number() == 1,elevated_value_mean,NA),
                normalized_expression = ifelse(row_number() == 1,normalized_expression,NA),
                log2_expression_value = ifelse(row_number() == 1,log2_expression_value,NA)) %>% 
  dplyr::ungroup()

## remake the information by sample
bio_rep_data <- raw_data2 %>% 
  dplyr::select(c(
    file_name,Pos,sample_name,zscoring_group,testing_group,
    Gene,type_of_gene,
    elevated_value_mean,normalized_expression,log2_expression_value,qc_flag)) %>% 
  dplyr::filter(!is.na(elevated_value_mean))

ttest_ready_data <- bio_rep_data
ttest_ready_data <- ttest_ready_data %>% 
  dplyr::filter(qc_flag == "GOOD")

bio_rep_data_red <- bio_rep_data %>% 
  dplyr::group_by(sample_name, Gene) %>% 
  dplyr:: mutate(normalized_expression_mean = mean(normalized_expression[qc_flag == "GOOD"]),
                 log2_expression_value_mean = mean(log2_expression_value[qc_flag == "GOOD"])) %>% 
  dplyr::mutate(normalized_expression_sd = sd(normalized_expression[qc_flag == "GOOD"]),
                log2_expression_value_sd = sd(log2_expression_value[qc_flag == "GOOD"])) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c(file_name,Pos,sample_name,zscoring_group,testing_group,
                  Gene, type_of_gene,
                  normalized_expression_mean, log2_expression_value_mean,
                  normalized_expression_sd, log2_expression_value_sd)) %>% 
  dplyr::distinct()

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
  dplyr::select(c(file_name,Pos,qc_flag))

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

raw_data2 <- raw_data2 %>% 
  dplyr::select(-qc_flag)

raw_data2 <- merge(
  x = raw_data2, y = bio_rep_data,
  by = c("file_name","Pos"), all.x = TRUE)
raw_data2 <- merge(
  x = raw_data2, y = bio_rep_data_red,
  by = c("file_name","Pos"), all.x = TRUE)

raw_data2 <- raw_data2 %>% 
  dplyr::arrange(replicate_name, Gene, is.na(elevated_value_mean))

openxlsx::write.xlsx(x = ttest_ready_data,
                     file = "./results/statistical_results2.xlsx")
openxlsx::write.xlsx(x = raw_data2, 
                     file = "./results/final_data2.xlsx")
raw_data <- raw_data2
remove(ttest_ready_data);remove(bio_rep_data);remove(bio_rep_data_red);remove(raw_data2)

