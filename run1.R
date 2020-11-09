# In this file we only set up the initial values for mu_t, mu_d, mu_w, mu_bg, f, g, mu0 and A

library(foreach)
library(doParallel)
library(polyCub)
library(spatstat) # must run R with option "--max-ppsize=100000"
gpclibPermit()

# functions
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

# SET UP KERNEL BW
# from paper
bw_daily <- 0.03
bw_weekly <- 0.5

# maximum trigger ranges allowed
max_t <- 15 # days
max_d <- 2 # km

# ------------------
# SETUP BASEVALUES
# ------------------

# 1. smoothing daily 

# daily.base is splitting the day into roughly 6min chunks
daily.base <- seq(0, 1, 0.005)

# calculate how far into the day each event happened. unit: 1/days
new.marks <- a$days - as.integer(a$days)

# get intensity at points throughout the day
temp <- hist.weighted(new.marks, rep(1, nrow(a)), breaks = daily.base)

# smooth points from midpoints density with bw 0.03
daily.basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_daily)

# standardize
daily.basevalue <- daily.basevalue/ mean(daily.basevalue)


# 2. smoothing weekly

# weekly.base is splitting the week into roughly 50min breaks
weekly.base <- seq(0, 7, 0.005)

# calculate how far into the week each event happened. unit: 1/days
new.marks <- a$days - as.integer(a$days / 7) * 7 

weights <- 1 / (as.integer(TT/7) + (a$days - as.integer(a$days / 7) * 7 > TT - as.integer(TT / 7) * 7)) 

# as before, get intensity at points throughout the week
temp <- hist.weighted(new.marks, weights, breaks = weekly.base)

# smooth density
weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_weekly)
# standardize
weekly.basevalue <- weekly.basevalue / mean(weekly.basevalue)


# 3. smoothing all trend

# split time interval into 3.65 days
trend.base <- seq(0, 730, 0.005)

trend.basevalue <- rep(0, length(time.marks))

# at each step, we add the density at the event, 
# normalized by density that lies between 0 and TT (our event window)
for(i in 1:nrow(a)){
  trend.basevalue <- (trend.basevalue + dnorm(a$days[i] - time.marks, 0, 100)
                      /(pnorm(TT, a$days[i], 100) - pnorm(0, a$days[i], 100)))    
}

# standardize
trend.basevalue <- trend.basevalue / mean(trend.basevalue)


# 4. Set up initial spatial background rate

# this background rates needs to integrate 1 over the study area

# define fine grid over entire study area
background.base <- list(x = seq(Xrange[1], Xrange[2], 0.002), 
                        y = seq(Yrange[1], Yrange[2], 0.002))

# create matrix of all locations on grid
background.basex <- background.base$x %o% rep(1, length(background.base$y))
background.basey <- rep(1, length(background.base$x))%o% background.base$y

# for each location, check if it's outside the boundary (-1), on the boundary (0) or inside (1)
background.marks <- matrix(inpoly(background.basex,   
                                  background.basey, 
                                  city.boundary$x, 
                                  city.boundary$y), 
                           ncol = length(background.base$y))

# predefine empty matrix
background.basevalue <- matrix(0, 
                               nrow = length(background.base$x), 
                               ncol = length(background.base$y))

if(!file.exists("Background.Smoothers")){
  system("mkdir Background.Smoothers")
}

# for each location on the grid, add density of each event
for(i in 1:nrow(a)){
  
  fn <- paste0("Background.Smoothers/", "bgsmoother-", i, ".val")
  
  if(!file.exists(fn)){
    
    bgsmoother <- (
      dnorm(background.basex, 
            a$coorx[i],
            a$bandwidth[i]) * 
      dnorm(background.basey,
            a$coory[i], 
            a$bandwidth[i]) /
      a$bg.integral[i]
    )
    save(bgsmoother,file=fn)
  } else{
    load(fn)
  }
  
  background.basevalue <- background.basevalue + bgsmoother
}

# standardize the background so its average inside the study window is 1
background.basevalue <- background.basevalue / mean(background.basevalue[background.marks >= 0])


# 5. Setting up triggering

# Temporal triggering

# these are just some initial values where we hardcode that we don't 
# expect any effect after more than 15 days
excite.temporal.basevalue <- (0.05 + seq(0, 15, 0.005)/10)^(-1.03)

# we normalize by the integral over the entire time window
excite.temporal.basevalue <- excite.temporal.basevalue / simpson(excite.temporal.basevalue, 0.005)



# spatial triggering

# again, we hardcode that there is no effect more than 2km away
excite.spatial.base.x <- seq(-2, 2, 0.002)
excite.spatial.base.y <- seq(-2, 2, 0.002)
d <- length(excite.spatial.base.x)

# not entirely sure what's happening here
excite.spatial.basevalue <- matrix(
  1/(abs(excite.spatial.base.x %o% rep(1, length(d)))^2 + 
       abs(rep(1, length(d)) %o% excite.spatial.base.y)^2 + 1), 
  ncol = d, 
  nrow = d)

# normalize by integral over the entire spatial window
excite.spatial.basevalue <- excite.spatial.basevalue / simpson.2D(excite.spatial.basevalue, 0.002, 0.002)


# Obtaining lambda and mu0

# now we can evaluate the lambda function 


# Define a couple of functions which (linearly) interpolate between x and y

trend.fun <- approxfun(time.marks, trend.basevalue, yleft=0, yright=0)

weekly.fun <- function(x){
  approxfun(weekly.base, weekly.basevalue,             
            yleft=0, yright=0)(x - as.integer(x/7)*7)
}

daily.fun <- function(x){
  approxfun(daily.base, daily.basevalue,
            yleft=0, yright=0)(x - as.integer(x))
}

# interpolate mu_b
background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x, 
                                                                 y=background.base$y, 
                                                                 z=background.basevalue),
                                                        loc=cbind(x=c(x), y=c(y))))

# initial guess for g(t)
excite.temporal.fun <- approxfun(seq(0, 15, 0.005)+0.6e-12, 
                                 excite.temporal.basevalue, 
                                 yleft=0, yright=0)

# interpolate spatial trigger component (initial guess for h(s))
excite.spatial.fun <- function(x,y){
  temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                  y=excite.spatial.base.y, 
                                  z=excite.spatial.basevalue),
                         loc=cbind(x=c(x), y=c(y))) 
  temp[is.na(temp)] <- 0
  temp
}

# evaluate background at all events
bgrates.at.events.no.mu <- (trend.fun(a$days) * weekly.fun(a$days) * daily.fun(a$days)
                            * background.spatial.fun(a$coorx,a$coory))

# mean of background over entire space
bgrates.at.all.no.mu <- (mean(trend.fun(time.marks) * weekly.fun(time.marks) * daily.fun(time.marks))* TT *
                         mean(background.spatial.fun(background.basex,background.basey) * background.marks)*
                         (Xrange[2] - Xrange[1])*(Yrange[2]-Yrange[1]))

# initalize triggered events
triggers.at.events.no.A <- rep(0, nrow(a))    
triggers.at.all.no.A <- 0

# I made some simplifications 
bg_weight <- sum(background.marks > 0)/length(background.marks)
K <- (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1])

# this loop is very slow
for (i in 1:nrow(a)){
  if (i %% 100 == 0) print(paste("on:", i))
  t_temp <- excite.temporal.fun(a$days - a$days[i]) * excite.spatial.fun(a$coorx - a$coorx[i], a$coory - a$coory[i])
  triggers.at.events.no.A <- triggers.at.events.no.A + t_temp
  
  temp <- excite.spatial.fun(background.basex - a$coorx[i], 
                             background.basey - a$coory[i])
  
  triggers.at.all.no.A <- triggers.at.all.no.A + mean(excite.temporal.fun(time.marks - a$days[i])) * TT * mean(temp[background.marks > 0]) * bg_weight * K
}


# Update A and mu for the first time
A <- 0.5
mu <- 0.7

# update according to equations (34) and (35)
NegLogLikehood <- function(x){
  mu <- x[1]^2
  A <- x[2]^2
  
  lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
  lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A
  
  - sum(log(lambda.at.events)) + lambda.at.all
}

res.optim <- optim(par=sqrt(c(A, mu)), NegLogLikehood, control=list(trace=6))

# update
mu <- res.optim$par[1]^2
A <- res.optim$par[2]^2


lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A

bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events

