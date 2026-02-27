expand_beta_blocks <- function(beta_values, beta_times, time_vec) {
  out <- numeric(length(time_vec))
  for (i in seq_along(beta_times)) {
    if (i < length(beta_times)) {
      idx <- which(time_vec >= beta_times[i] & time_vec < beta_times[i+1])
    } else {
      idx <- which(time_vec >= beta_times[i])
    }
    out[idx] <- beta_values[i]
  }
  out
}
