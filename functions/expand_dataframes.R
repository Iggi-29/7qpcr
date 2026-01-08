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
