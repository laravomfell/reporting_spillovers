library(spatstat)
library(sf)
library(polyCub)
library(mvtnorm)
gpclibPermit()

library(foreach)
library(doParallel)

# set seed to jitter
set.seed(1)

no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# pretty self explanatory
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}

# read crime data
da <- read.csv(file="da.csv")

da$datetime_unif <- as.POSIXct(da$datetime_unif, tz = "Europe/London")

# calculate time in seconds since t0, turn into 1/day units
da$days <- (unclass(da$datetime_unif) - unclass(min(as.POSIXct(as.Date(da$datetime_unif)))))/24/60/60

da <- da[order(da$days),]

# turn coordinates into km
da$coorx <- da$X / 1000
da$coory <- da$Y / 1000

# add jitter
da$coorx <- da$coorx + runif(nrow(da), -0.0005, 0.0005)
da$coory <- da$coory + runif(nrow(da), -0.0005, 0.0005)

# set bandwidth around events
da$bandwidth <- rep(0.01, nrow(da))

for(i in 1:nrow(da)){
   # calculate sq distances to all points not == i
   temp <- dist.squared(da$coorx[i], da$coory[i], da$coorx[-i], da$coory[-i])
   # n_p = 5
   temp2 <- sqrt(sort(temp[temp >= 0.002^2])[5])
   # replace bandwidth if too small otherwise
   if(da$bandwidth[i] <= temp2){
      da$bandwidth[i] <- temp2
   }
}

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


foo <- function(x, mu, sigma){
   mvtnorm::dmvnorm(x, mean = mu, sigma = sigma * diag(length(mu)), checkSymmetry = FALSE)
}

da$bg_integral <- foreach(i = 1:nrow(da),
                          .combine = "c") %dopar% polyCub::polyCub.SV(w,
                                                                      foo,
                                                                      mu=c(da$coorx[i], da$coory[i]),
                                                                      sigma=sqrt(da$bandwidth[i]),
                                                                      plot=F,
                                                                      nGQ = 15)

write.csv(da, file = "da_type.csv", row.names = F)
stopCluster(cl)
