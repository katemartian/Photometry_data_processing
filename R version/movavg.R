movavg <- function(v, win)
{
  # Calculates moving average
  # 
  # Input
  #     v: vector to smooth
  #     win: window for moving average
  # Output
  #     smoothed signal
  #
  
  f = rep(1, win) / win
  n = length(v)
  v1 = c(rep(v[1], win), v, rep(v[n],win))
  v1 = stats::filter(v1, f)
  return(v1[win+1:n])
}
