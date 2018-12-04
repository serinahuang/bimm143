add <- function(x, y = 1) {
  # sum inputs of x and y
  x + y
}

add( c(1,2,3) )
add( c(1,2,3), 4)
add(1,2,2) # error
add(x=1, y="b") # error

# making a function
rescale <- function(x){
  rng <- range(x)
  (x - rng[1]) / (rng[2] - rng[1])
}

rescale(c(3:11))
rescale(c(8:10))
rescale(c(2,7))
rescale( c(1,2,NA,3,10) )
rescale( c(1,10,"string") )

# fix NA
range( c(1,2,NA,4), na.rm = TRUE)

rescale2 <- function(x, plot = TRUE){
  rng <- range(x, na.rm = TRUE)
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  if(plot){
    plot(answer, type = "b")
  }
}

# or this:
# rescale2 <- function(x, na.rm = TRUE){
# rng <- range(x, na.rm = na.rm)
# (x - rng[1]) / (rng[2] - rng[1])
# }

rescale( c(1,2,NA,3,10) )
rescale2( c(1,2,NA,3,10) ) # shows missing dot for NA

rescale3 <- function(x, na.rm = TRUE, plot = FALSE){
  if(na.rm){
    rng <- range(x, na.rm = TRUE)
  }
  else{
    rng <- range(x) # do you need this line?
  }
  print("Hello, here is your range:")
  print(rng)
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  print("Here is your answer:")
  print(answer)
  print("Is it me you're looking for?")
  
  if(plot){
    plot(answer, type = "o", lwd = 4, col = rainbow(7), pch = 15)
  }
}

rescale3( c(1,2,4,5,6,NA,8), plot = TRUE )

## 1B. Improve code for analyzing protein drug interactions
library(bio3d)

x = "4AKE"

# the drug_effect function
drug_effect <- function(x){
  s <- read.pdb(x)
  chainA <- trim.pdb(s, chain="A", elety="CA")
  b <- chainA$atom$b
  plotb3(b, sse=chainA, typ="l", ylab="B factor", col="black")
  points(b, typ="l", col="red") # overlay
}

# to turn off black and gray rectangles, which represent secondary structure
plotb3(s3.b, sse=NULL, typ="l", ylab="Bfactor")

drug_effect(x)
