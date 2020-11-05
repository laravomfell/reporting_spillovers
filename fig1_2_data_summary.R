library(foreach)
library(doParallel)
library(polyCub)
library(spatstat) # must run R with option "--max-ppsize=100000"
gpclibPermit()

# functions
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}
source("utils.R")

# read data (from source("prepare_data.R"))
a <- read.table("type1_crime_data.table")
# BOUNDARY of CASTELLON
city.boundary <- read.csv(file = "castallon_city_boundary.csv")
city.boundary <- list(x = city.boundary$X, y = city.boundary$Y)

city.boundary$x <-  city.boundary$x /1000
city.boundary$y <- city.boundary$y /1000

city.boundary$x <- rev(city.boundary$x)
city.boundary$y <- rev(city.boundary$y)

# DEFINE PP PARS
# TT = length of the time interval
# TT <- 730
TT <- as.numeric(as.Date(max(a$date)) - as.Date(min(a$date))) + 1

# time stamps as 1 year/100 (= 3.65 days)
time.marks <- seq(0, TT, 1/(TT/3.65))

# range of coordinates (not quite, a little manual)
Xrange <- c(749.440, 754.264)
Yrange <- c(4428.644, 4432.570)

# a contains the crime information
# city.boundary is a list of X and Y boundary coords
par(mfrow=c(2,2),cex=1.5)

plot(a$coorx, a$coory,col= rainbow(80)[a$days/10],xlab="X", ylab="Y",cex=0.3, main="(a)")
points(city.boundary$x, city.boundary$y, type="l", lwd=2, col=gray(0.5))

plot(a$days, a$coory, xlab="Days", col= rainbow(80)[a$days/10],ylab="Y",cex=0.3, main="(b)")
plot(a$coorx, a$days, xlab="X", col= rainbow(80)[a$days/10],ylab="Days",cex=0.3, main="(c)")

plot(c(0,a$days,TT), c(0:(nrow(a)+1)) , type="s", lwd=2, xlab="Days", ylab="Cum. Freq.", main="(d)")