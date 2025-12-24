##########################################
### Work qPCR data in an automated way ###
##########################################

### WD
setwd("/home/labs/eclab/ijarne/Desktop/")

### libraries
library(readr)
library(dplyr)
library(purrr)

### Constants
raw_data_file <- "./7qpcr/raw_data_to_play/20251212_SPOCK1_DN-TRES.txt"
samples_locations <- ""
genes_locations <- ""

### Data importation
raw_data <- readr::read_tsv(
  file = raw_data_file, 
  skip = 1, 
  locale = readr::locale(decimal_mark = ","))

### Genes work
## Genes location
# Based on the row location - in the future we could pattern recognition maybe
genes_location <- list("RHOA" = c("A","B","I","J"),
                       "SPOCK" = c("C","D","K","L"),
                       "AXL" = c("E","F","M","N"),
                       "PMEL" = c("G","H","O","P"))

genes_location_df <- data.frame(
  gene = rep(names(genes_location), lengths(genes_location)),
  row_name  = unlist(genes_location),
  stringsAsFactors = FALSE
)
remove(genes_location)

### Housekeeping genes
genes_housekeeping <- "RHOA"


### Samples work
# samples_location <- ""
samples_location_74 <- 
  list(
    ### NORMAL
    "MM074.DN.1" = 
      list(
        cols = c("1","2","3"),
        rows = c("A","C","E","G")),
    "MM074.DN.2" = 
      list(
        cols = c("1","2","3"),
        rows = c("J","L","N","P")),
    "MM074.DN.3" = 
      list(
        cols = c("4","5","6"),
        rows = c("A","C","E","G")),
    "MM074.DN.4" = 
      list(
        cols = c("7","8","9"),
        rows = c("A","C","E","G")),
    "MM074.DN.5" = 
      list(
        cols = c("10","11","12"),
        rows = c("A","C","E","G")),
    ### RES
    "MM074.R.1" = 
      list(
        cols = c("13","14","15"),
        rows = c("A","C","E","G")),
    "MM074.R.2" = 
      list(
        cols = c("16","17","18"),
        rows = c("A","C","E","G")),
    "MM074.R.3" = 
      list(
        cols = c("19","20","21"),
        rows = c("A","C","E","G")),
    "MM074.R.4" = 
      list(
        cols = c("22","23","24"),
        rows = c("A","C","E","G")))


samples_location_99 <- 
  list(
    ### NORMAL
    "MM099.DN.1" = 
      list(
        cols = c("1","2","3"),
        rows = c("I","K","M","O")),
    "MM099.DN.2" = 
        list(
          cols = c("4","5","6"),
          rows = c("I","K","M","O")),
    "MM099.DN.3" = 
      list(
        cols = c("4","5","6"),
        rows = c("J","L","N","P")),
    "MM099.DN.4" = 
      list(
        cols = c("7","8","9"),
        rows = c("I","K","M","O")),
    "MM099.DN.5" = 
      list(
        cols = c("7","8","9"),
        rows = c("J","L","N","P")),
    "MM099.DN.6" = 
      list(
        cols = c("10","11","12"),
        rows = c("I","K","M","O")),
    ### RES
    "MM099.R.1" = 
      list(
        cols = c("13","14","15"),
        rows = c("I","K","M","O")),
    "MM099.R.2" = 
      list(
        cols = c("16","17","18"),
        rows = c("I","K","M","O")),
    "MM099.R.3" = 
      list(
        cols = c("16","17","18"),
        rows = c("J","L","N","P")),
    "MM099.R.4" = 
      list(
        cols = c("19","20","21"),
        rows = c("I","K","M","O")),
    "MM099.R.5" = 
      list(
        cols = c("19","20","21"),
        rows = c("J","L","N","P")),
    "MM099.R.6" = 
      list(
        cols = c("22","23","24"),
        rows = c("I","K","M","O")))
        

samples_location_383 <- 
  list(
    ### NORMAL
    "MM0383.DN.1" = 
      list(
        cols = c("1","2","3"),
        rows = c("B","D","F","H")),
    "MM0383.DN.2" = 
      list(
        cols = c("4","5","6"),
        rows = c("B","D","F","H")),
    "MM0383.DN.3" = 
      list(
        cols = c("7","8","9"),
        rows = c("B","D","F","H")),
    "MM0383.DN.4" = 
      list(
        cols = c("10","11","12"),
        rows = c("B","D","F","H")),
    "MM0383.DN.5" = 
      list(
        cols = c("10","11","12"),
        rows = c("J","L","N","P")),
    "MM0383.DN.6" = 
      list(
        cols = c("13","14","15"),
        rows = c("J","L","N","P")),
    ### RES
    "MM0383.R.1" = 
      list(
        cols = c("13","14","15"),
        rows = c("B","D","F","H")),
    "MM0383.R.2" = 
      list(
        cols = c("16","17","18"),
        rows = c("B","D","F","H")),
    "MM0383.R.3" = 
      list(
        cols = c("19","20","21"),
        rows = c("B","D","F","H")),
    "MM0383.R.4" = 
      list(
        cols = c("22","23","24"),
        rows = c("B","D","F","H"))
    
    )

### final list
sample_locations <- c(samples_location_74, samples_location_99, samples_location_383)

### Generate a df with sample locations
sample_locations_df <- do.call(
  rbind,
  lapply(names(sample_locations), function(s) {
    x <- sample_locations[[s]]
    data.frame(
      sample   = s,
      Pos = as.vector(outer(x$rows, x$cols, paste0)),
      stringsAsFactors = FALSE
    )
  })
)
remove(samples_location_383);remove(samples_location_74);remove(samples_location_99)
remove(sample_locations)

#### mutate raw_data
raw_data <- raw_data %>%
  # get the useful data
  dplyr::select(c(Pos,Cp)) %>%
  # generate col and row data
  dplyr::rowwise() %>% 
  dplyr::mutate(row_name = strsplit(x = Pos, split = "")[[1]][1],
                col_name = strsplit(x = Pos, split = "")[[1]][2]) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(row_name, col_name, .before = Pos) %>% 
  # Add the gene info
  dplyr::mutate(Gene = genes_location_df$gene[match(row_name, genes_location_df$row_name)]) %>% 
  dplyr::mutate(type_of_gene = ifelse(Gene %in% genes_housekeeping,"HK","Normal")) %>% 
  # Add replicate_name
  dplyr::mutate(replicate_name = sample_locations_df$sample[match(Pos, sample_locations_df$Pos)])

remove(genes_location_df);remove(sample_locations_df);remove(genes_housekeeping);remove(raw_data_files)

raw_data <- raw_data %>% 
  dplyr::select(-c(row_name, col_name)) %>% 
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

##### 2 remove samples with HK35
raw_data <- raw_data %>% 
  dplyr::group_by(replicate_name) %>%
  dplyr::mutate(
    good_replicate = ifelse(
      any(type_of_gene == "HK" & Cp >= 35, na.rm = TRUE),
      "bad","good")) %>% 
  dplyr::ungroup()

##### 3 NAs to 35
raw_data <- raw_data %>%
  dplyr::mutate(Cp = ifelse(is.na(Cp),35,Cp))

##### 4 check for bad values
raw_data <- raw_data %>%
  group_by(replicate_name, Gene) %>%
  mutate(
    comparison_results = map_chr(Cp, ~ {
      comps <- abs(.x - Cp) <= 0.5
      comps <- comps[!is.na(comps)]
      comps <- comps[-which(Cp == .x)[1]]  # remove self-comparison

      paste0(ifelse(comps, "good", "bad"), collapse = ",")
    })
  ) %>%
  ungroup() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(comparison_results2 = paste0(x = sort(unlist(strsplit(comparison_results, split = ",")[[1]])), collapse = ",")) %>% 
  dplyr::ungroup()

##### 5 delta Ct calculatio
# raw_data <- raw_data %>%
#     dplyr::group_by(replicate_name, Gene) %>%
#     dplyr::mutate(
#       n_good = sum(good_value == "good", na.rm = TRUE),
#       mean_value = dplyr::case_when(
#         n_good >= 2 ~ mean(Cp[good_value == "good"], na.rm = TRUE),
#         n_good <= 1 ~ mean(Cp, na.rm = TRUE)),
#       sd_value = dplyr::case_when(
#         n_good >= 2 ~ sd(Cp[good_value == "good"], na.rm = TRUE),
#         n_good <= 1 ~ sd(Cp, na.rm = TRUE),
#         TRUE        ~ NA_real_
#         )) %>%
#     dplyr::ungroup()
  

##### 6 delta delta CT value
raw_data <- raw_data %>%
  group_by(replicate_name) %>%
  mutate(
    hk_mean_value = mean(mean_value[type_of_gene == "HK"],
                         na.rm = TRUE),
    deltaCT = mean_value/hk_mean_value) %>%
  dplyr::ungroup()

##### Delta delta CT
### 5 0.5^Cp
raw_data <- raw_data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(elevated_value = 0.5^Cp) %>%
  dplyr::ungroup()

##### 6 mean of the CT
raw_data <- raw_data %>%
  dplyr::group_by(replicate_name, Gene) %>%
  dplyr::mutate(
    n_good = sum(good_value == "good", na.rm = TRUE),
    mean_value = dplyr::case_when(
      n_good >= 2 ~ mean(elevated_value[good_value == "good"], na.rm = TRUE),
      n_good <= 1 ~ mean(elevated_value, na.rm = TRUE)),
    sd_value = dplyr::case_when(
      n_good >= 2 ~ sd(elevated_value[good_value == "good"], na.rm = TRUE),
      n_good <= 1 ~ sd(elevated_value, na.rm = TRUE),
      TRUE        ~ NA_real_
      )) %>%
  dplyr::ungroup()


##### 7 delta CT
raw_data <- raw_data %>%
  group_by(replicate_name) %>%
  mutate(
    hk_mean_value = mean(mean_value[type_of_gene == "HK"],
                         na.rm = TRUE),
    deltaCT = mean_value/hk_mean_value) %>%
  dplyr::ungroup()

##### 8 Clean for Excel
raw_data <- raw_data %>%
  group_by(replicate_name, Gene) %>%
  mutate(
    row_in_group = row_number(),
    mean_value   = ifelse(row_in_group == 1, mean_value, NA_real_),
    sd_value     = ifelse(row_in_group == 1, sd_value, NA_real_),
    deltaCT      = ifelse(row_in_group == 1, deltaCT, NA_real_)
  ) %>%
  ungroup() %>%
  select(-row_in_group) %>%
  dplyr::arrange(replicate_name,type_of_gene,Gene)

openxlsx::write.xlsx(file = "./final_excel_for_L.xlsx", x = raw_data)



