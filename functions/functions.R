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
# run test function
run_test <- function(df, value_col, group_col, test_type = "ttest") {
  x <- df[[value_col]]
  g <- df[[group_col]]
  
  if (length(unique(g)) < 2 || length(x) < 2)
    return(NA_real_)
  
  if (test_type == "ttest") {
    t.test(x ~ g)$p.value
  } else {
    wilcox.test(x ~ g)$p.value
  }
}
## detect frontiers
detect_frontiers <- function(x, magnitude_jump = 2) {
  x <- as.numeric(x)
  
  # calculate the order of magnitude of each value
  orders <- floor(log10(x))
  
  # difference between consecutive orders
  diff_orders <- abs(diff(orders))
  
  # detect jumps >= magnitude_jump
  frontier_indices <- which(diff_orders >= magnitude_jump) + 1  # +1 because diff shifts by 1
  
  return(frontier_indices)
}
