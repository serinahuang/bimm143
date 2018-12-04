map.colors2 <- function (x, high.low, palette) {
  # determine percent values of the 'high.low' range
  percent <- (x - high.low[1]) / (high.low[2] - high.low[1])
 
  # find corresponding index position in the color 'palette'
  # +1 is if length = 0
  index <- round ( length(palette) * percent ) + 1
  
  return (palette[index])
}


map.colors3 <- function (x,
                         low.high = range(x),
                         palette = cm.colors(100)) {
  # determine percent values of the 'high.low' range
  percent <- (x - high.low[1]) / (high.low[2] - high.low[1])
  
  # find corresponding index position in the color 'palette'
  # +1 is if length = 0
  index <- round ( length(palette) * percent ) + 1
  
  return (palette[index])
}
