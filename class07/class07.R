## Revisit functions

source("http://tinyurl.com/rescale-R")

# Error messages to user are helpful
# Come up with application for each example

## Get working snippet
# Start with a simple steup to which we know the answer
x <- c(1, 2, NA, 4, NA)
y <- c(NA, 3, NA, 3, 4)

is.na(x)
is.na(y)
is.na(x) & is.na(y)

# True and false are binary values, i.e. they can be summed!
sum(is.na(x))
sum(is.na(y))
sum( is.na(x) & is.na(y) )

## Turn working snippet into a function
both_na <- function(x,y) {
  sum( is.na(x) & is.na(y) )
}

## Test the function by breaking it in more complicated cases
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)
both_na(x, y1)
both_na(x, y2)
# why does it do this?
# recall recycling! Shorter vector will recycle its value to match the longer vector
x2 <- c(1, NA, NA)
both_na(x2, y2)
# if you still don't believe, see with rbind
rbind(x2, y2)
rbind(x, y2)
rbind(x, c(y2, y2))

# but this is the wrong result! Fix function to include stop message
both_na2 <- function(x,y) {
  
  if (length(x) != length(y)) {
    stop("Input x and y should be of same length!")
  }
  
  sum( is.na(x) & is.na(y) )
}

both_na2(x, y2)

c(F, F, T, F)
which( c(F, F, T, F) ) # tells which element is true

res <- both_na3(x, y1)
res$which
res$number

# for commenting, document WHY you're doing it, not WHAT


## Intersecting!
df1
df2

x <- df1$IDs
y <- df2$IDs
x
y
intersect(x,y)
# read help page and found See Also
inds <- x %in% y
x[inds]
# WOW!
y %in% x
y [y %in% x]

df1.inds <- df1$IDs %in% df2$IDs
exp_1 <- df1[df1.inds, ]

df2.inds <- df2$IDs %in% df1$IDs
exp_2 <- df2[df2.inds, 2]
cbind(exp_1, exp_2)
