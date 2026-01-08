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
