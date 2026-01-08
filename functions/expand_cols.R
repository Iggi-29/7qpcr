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
