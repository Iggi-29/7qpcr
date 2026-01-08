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
