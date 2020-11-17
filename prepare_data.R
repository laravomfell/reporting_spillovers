library(spatstat)
library(sf)
library(polyCub)
library(mvtnorm)
gpclibPermit()

library(foreach)
library(doParallel)

no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# pretty self explanatory
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}

# read crime data
da <- read.csv(file="da_small.csv")
# keep only the original reports
da <- da[da$ms_type == 0,]

da$datetime_unif <- as.POSIXct(da$datetime_unif, tz = "Europe/London")

# calculate time in seconds since t0, turn into 1/day units
da$days <- (unclass(da$datetime_unif) - unclass(min(da$datetime_unif)))/24/60/60

da <- da[order(da$days),]

# turn coordinates into km
da$coorx <- da$X / 1000
da$coory <- da$Y / 1000

# don't add jitter
# a$coorx =a$coorx /1000+runif(nrow(a), -0.0005, 0.0005)
# a$coory = a$coory /1000+runif(nrow(a), -0.0005, 0.0005)

# set bandwidth around events
da$bandwidth <- rep(0.01, nrow(da))

for(i in 1:nrow(da)){
   # calculate sq distances to all points not == i
   temp <- dist.squared(da$coorx[i], da$coory[i], da$coorx[-i], da$coory[-i])
   # keep fifth smallest distance above threshold???
   # maybe something like we need at least five points nearby??
   temp2 <- sqrt(sort(temp[temp >= 0.002^2])[5])
   # replace bandwidth if too small otherwise
   if(da$bandwidth[i] <= temp2){
      da$bandwidth[i] <- temp2
   }
}


# Read boundary
boundary <- read_sf("cov.shp", crs = 27700)
boundary <- st_union(boundary)
boundary <- st_as_sf(boundary)
# fix tiny hole
boundary <- st_buffer(boundary, 1)
# simplify boundary after union
boundary <- st_simplify(boundary, preserveTopology = TRUE, dTolerance = 0.05)
# extract boundaries
bbox <- st_bbox(boundary) / 1000

boundary <- data.frame(st_coordinates(boundary)[, c("X", "Y")])
boundary$x <- boundary$X / 1000
boundary$y <- boundary$Y / 1000

# reverse x and y for owin
boundary$x <- rev(boundary$x)
boundary$y <- rev(boundary$y)

da$bg_integral <- rep(0, nrow(da))

w <- owin(c(bbox["xmin"], bbox["xmax"]), c(bbox["ymin"], bbox["ymax"]), poly = boundary)

#for(i in 1:nrow(da)){
#   if (i %% 100 == 0) print(paste("on:", i))
#   # calculate exact integral
#   da$bg_integral[i] <- polyCub.exact.Gauss(w,
#                                           mean=c(da$coorx[i], da$coory[i]),
#                                           Sigma=diag(da$bandwidth[i], 2), 
#                                           plot=F)
# }

write.table(a, file="type1_crime_data.table")