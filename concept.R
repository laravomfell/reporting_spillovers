# This slightly modified code adapted from Zhuang/Mateu proceeds as follows:

# SETUP
# initialize mu_daily, mu_weekly, mu_trend, mu_bg
# initalize g(t) and h(s)

# get background = mu_daily * mu_weekly * mu_trend * mu_bg
# we want both background_i (at the events)
# and the integral over background
# \int_0^T \int_S mu_daily * mu_weekly * mu_trend * mu_bg dt ds

# trigger = g(t) * h(s)
# again, both as trigger_i (at events)
# and integral over trigger
# sum_i \int_t_i ^T \int_S g(t-t_i)h(s-s_i) dt ds

# then update A and mu_0 based on loglik

# then calculate lambda_i = mu0 * background_i + A * trigger_i
# and lambda = mu0 * background_all + A * trigger_all

# now all quantities are initialized.

# to later update g(t) and h(s) we need two quantities which we can precalculate now.

# 1. repetition terms calculating how often the trigger function is potentially "triggered"/ is applicable
# calculate \sum_i I(t_i + t \leq T) and \sum_i I(s_i + s \in s)

# 2. a precursor of the rho_ij matrix giving us "candidate" pairs of (i,j) which are close enough in space/time
# to be potentially triggered 

# this concludes the setup.

# ENTER LOOP:

# get weights wi_daily, wi_weekly, wi_trend, wi_bg

# use those to get mu_daily, mu_weekly, mu_trend, mu_bg

# get g(t) propto equation:
# calculate g(t) edge correction
# get g(t)

# get h(s) propto part
# calculate h(s) edge correction
# get h(s)

# get background_i and int_background

# get trigger_i and int_trigger

# then update A and mu_0

# then calculate lambda_i = mu0 * background_i + A * trigger_i
# and lambda = mu0 * background_all + A * trigger_all

# terminate or loop back

library(spatstat) # must run R with option "--max-ppsize=100000"
library(purrr)

# functions
source("utils.R")

# read data (from source("prepare_data.R"))
a <- read.table("type1_crime_data.table")
# BOUNDARY of CASTELLON
city.boundary <- read.csv(file = "castallon_city_boundary.csv")
city.boundary <- list(x = city.boundary$X, y = city.boundary$Y)

city.boundary$x <- city.boundary$x /1000
city.boundary$y <- city.boundary$y /1000

city.boundary$x <- rev(city.boundary$x)
city.boundary$y <- rev(city.boundary$y)

# DEFINE PP PARS
# TT = length of the time interval
# TT <- 730
TT <- as.numeric(as.Date(max(a$date)) - as.Date(min(a$date))) + 1

# time stamps as 1 year/100 (= 3.65 days)
time_marks <- seq(0, TT, 1/(TT/3.65))

# range of coordinates (not quite, a little manual)
x_range <- c(749.440, 754.264)
y_range <- c(4428.644, 4432.570)

# range of area
ra <- (x_range[2]-x_range[1])*(y_range[2]-y_range[1])

# SET UP KERNEL BW
# from paper
bw_daily <- 0.03
bw_weekly <- 0.5

# maximum trigger ranges allowed
max_t <- 15 # days
max_s <- 2 # km

time_units <- 0.005
space_units <- 0.002

# STEP 1

# init mu_d
# 1. smoothing daily 

# daily.base is splitting the day into roughly 6min chunks
daily_base <- seq(0, 1, time_units)

# calculate how far into the day each event happened. unit: 1/days
day_marks <- a$days - as.integer(a$days)

# get intensity at points throughout the day
temp <- hist.weighted(day_marks, rep(1, nrow(a)), breaks = daily_base)

# smooth points from midpoints density with bw 0.03
daily_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_daily)

# standardize
daily_basevalue <- daily_basevalue/ mean(daily_basevalue)

# init mu_w
# weekly_base is splitting the week into roughly 50min breaks
weekly_base <- seq(0, 7, time_units)

# calculate how far into the week each event happened. unit: 1/days
week_marks <- a$days - as.integer(a$days / 7) * 7 

week_weights <- 1 / (as.integer(TT/7) + (a$days - as.integer(a$days / 7) * 7 > TT - as.integer(TT / 7) * 7)) 

# as before, get intensity at points throughout the week
temp <- hist.weighted(week_marks, week_weights, breaks = weekly_base)

# smooth density
weekly_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_weekly)
# standardize
weekly_basevalue <- weekly_basevalue / mean(weekly_basevalue)

# init mu_t
# 3. smoothing all trend

# split time interval into 3.65 days
trend_base <- seq(0, 730, time_units)

# at each step, we add the density at the event, 
# normalized by density that lies between 0 and TT (our event window)
trend_basevalue <- map(a$days, function(x) dnorm(x - time_marks, 0, 100) / (pnorm(TT, x, 100) - pnorm(0, x, 100)))
trend_basevalue <- reduce(trend_basevalue, `+`)

# standardize
trend_basevalue <- trend_basevalue / mean(trend_basevalue)

# init mu_bg
# this background rates needs to integrate 1 over the study area

# define fine grid over entire study area
background_base <- list(x = seq(x_range[1], x_range[2], space_units), 
                        y = seq(y_range[1], y_range[2], space_units))

# create matrix of all locations on grid
background_basex <- background_base$x %o% rep(1, length(background_base$y))
background_basey <- rep(1, length(background_base$x))%o% background_base$y

# for each location, check if it's outside the boundary (-1), on the boundary (0) or inside (1)
background_marks <- matrix(inpoly(background_basex,   
                                  background_basey, 
                                  city.boundary$x, 
                                  city.boundary$y), 
                           ncol = length(background_base$y))

# predefine empty matrix
background_basevalue <- matrix(0, 
                               nrow = length(background_base$x), 
                               ncol = length(background_base$y))

if(!file.exists("Background.Smoothers")){
  system("mkdir Background.Smoothers")
}

# for each location on the grid, add density of each event
for(i in 1:nrow(a)){
  
  fn <- paste0("Background.Smoothers/", "bgsmoother-", i, ".val")
  
  if(!file.exists(fn)){
    
    bgsmoother <- (
      dnorm(background_basex, 
            a$coorx[i],
            a$bandwidth[i]) * 
        dnorm(background_basey,
              a$coory[i], 
              a$bandwidth[i]) /
        a$bg.integral[i]
    )
    save(bgsmoother,file=fn)
  } else{
    load(fn)
  }
  
  background_basevalue <- background_basevalue + bgsmoother
}

# standardize the background so its average inside the study window is 1
background_basevalue <- background_basevalue / mean(background_basevalue[background_marks >= 0])



# init g(t)
# Temporal triggering

# these are just some initial values where we don't 
# expect any effect after more than 15 days
g_base <- seq(0, max_t, time_units)
g_basevalue <- (0.05 + g_base/10)^(-1.03)

# we normalize by the integral over the entire time window
g_basevalue <- g_basevalue / simpson(g_basevalue, time_units)

# init h(s)
# spatial triggering

# there is no effect more than 2km away
h_base_x <- seq(-max_s, max_s, space_units)
h_base_y <- seq(-max_s, max_s, space_units)
d <- length(h_base_x)

# not entirely sure what's happening here
h_basevalue <- matrix(
  1/(abs(h_base_x %o% rep(1, d))^2 + 
       abs(rep(1, d) %o% h_base_y)^2 + 1), 
  ncol = d, 
  nrow = d)

# normalize by integral over the entire spatial window
h_basevalue <- h_basevalue / simpson.2D(h_basevalue, space_units, space_units)


# init mu0 & A
source("interpolate_fun.R")

# evaluate background at all events
bg_at_events_no_mu <- (trend_fun(a$days) * 
                       weekly_fun(a$days) * 
                       daily_fun(a$days) *
                       background_fun(a$coorx, a$coory))

# mean of background over entire space
bg_at_all_no_mu <- (mean(trend_fun(time_marks) * 
                         weekly_fun(time_marks) * 
                         daily_fun(time_marks)) * 
                    TT *
                    mean(background_fun(background_basex,
                                        background_basey) * 
                         background_marks) *
                    ra)

# trigger prob rho
rho_at_events_no_A <- map(1:nrow(a), 
                          function(i) g_fun(a$days - a$days[i]) * 
                                      h_fun(a$coorx - a$coorx[i], 
                                            a$coory - a$coory[i]))
rho_at_events_no_A <- reduce(rho_at_events_no_A, `+`)

# this loop takes approx 3h because we need to evaluate 
# h_fun repeatedly over the entire grid
bg_weight <- sum(background_marks > 0)/length(background_marks)

# get mean effect of g_fun over the entire time
# multiply through with all constants
rho_at_all_no_A <- map(a$days, function(x) mean(g_fun(time_marks - x)) * TT * ra)

for (i in 1:nrow(a)){
  if (i %% 100 == 0) print(paste("on:", i))
  temp <- h_fun(background_basex - a$coorx[i], 
                background_basey - a$coory[i])
  rho_at_all_no_A[[i]] <- rho_at_all_no_A[[i]] * mean(temp[background_marks > 0]) * bg_weight
}

rho_at_all_no_A <- reduce(rho_at_all_no_A, sum)


# Update A and mu for the first time
A <- 0.03
mu0 <- 0.77

# update according to equations (34) and (35)
neg_loglik <- function(x){
  mu0 <- x[1]^2
  A <- x[2]^2
  
  lambda_at_events <- mu0 * bg_at_events_no_mu + A * rho_at_events_no_A
  lambda_at_all <- mu0 * bg_at_all_no_mu + A * rho_at_all_no_A
  
  - sum(log(lambda_at_events)) + lambda_at_all
}

res <- optim(par=sqrt(c(mu0, A)), neg_loglik, control=list(trace=6))

# update mu0 and A
mu0 <- res$par[1]^2
A <- res$par[2]^2


# update other quantities
lambda_at_events <- mu0 * bg_at_events_no_mu + A * rho_at_events_no_A
lambda_at_all <- mu0 * bg_at_all_no_mu + A * rho_at_all_no_A

# calculate bg_probs (phi_i)
bg_probs <- mu0 * bg_at_events_no_mu / lambda_at_events


# pre-calculate the trigger terms

# candidate (i,j) pairs

# create a matrix ij.mat where 
# [, i] is the i index
# [, j] lists all events j where j < i
# and (i,j) are no more than 2km and 15 days apart
temp.mat <- (1:nrow(a)) %o% rep(1, nrow(a))

ij.mat <- cbind(c(t(temp.mat)), c(temp.mat))

ij.mat <- ij.mat[a$days[ij.mat[,1]] > a$days[ij.mat[,2]] & 
                   a$days[ij.mat[,1]] <= a$days[ij.mat[,2]] + 15.0 & 
                   abs(a$coorx[ij.mat[,1]] - a$coorx[ij.mat[,2]]) <=2 & 
                   abs(a$coory[ij.mat[,1]] -a$coory[ij.mat[,2]]) <=2,]


# repetance = "repetition corrections, i.e. for how many times the triggering effect at time lag t or the
# spatial offset (x, y) is observed"
# expand h_base
h_base_x <- h_base_x %o% rep(1, d)
h_base_y <- rep(1, d) %o% h_base_y

g_rep <- rep(0, length(g_base))
h_rep <- matrix(0, 
                ncol = length(h_base_x),
                nrow = length(h_base_y)) 


for(i in 1:nrow(a)){
  g_rep[g_base < TT - a$days[i]] <- g_rep[g_base < TT -a$days[i]] + 1  
  
  if(!file.exists("Excite.Spatail.Marks")){system("mkdir Excite.Spatail.Marks")}
  
  fn <- paste0("Excite.Spatail.Marks/crime1-", substr(10000 + i, 2, 5), ".mark")
  
  if(!file.exists(fn)){
    h_mark_temp <- matrix((inpoly(h_base_x + a$coorx[i], 
                                  h_base_y + a$coory[i], 
                                  city.boundary$x, 
                                  city.boundary$y) >=0),
                          ncol=length(h_base_x))
    save(h_mark_temp, file=fn)
  } else {
    load(fn)
  }
  
  h_rep <- h_rep + h_mark_temp
}

g_rep_fun <- approxfun(g_base, g_rep, yleft=1, yright=1)

h_rep_fun <- function(x,y){
  temp <- interp.surface(obj=list(x = h_base_x, 
                                  y = h_base_y, 
                                  z = h_rep), 
                         loc=cbind(x=c(x), y=c(y)))
  temp[is.na(temp)] <- 0
  temp
}


# ENTER THE LOOP ###############################################################


# get weights wi_daily, wi_weekly, wi_trend, wi_bg
# calculate w_i_t
wi_trend <- trend_fun(a$days)/lambda_at_events

# calculate w_i_d
wi_daily <- daily_fun(a$days) * background_fun(a$coorx, a$coory) / lambda_at_events

# calculate w_i_w
wi_weekly <- weekly_fun(a$days) * background_fun(a$coorx, a$coory) / lambda_at_events

# use those to get mu_daily, mu_weekly, mu_trend, mu_bg

# ESTIMATE MU_TREND

# same as before, just now including wi_trend
trend_basevalue <- map(a$days, function(x) wi_trend[i] * dnorm(x - time_marks, 0, 100) / (pnorm(TT, x, 100) - pnorm(0, x, 100)))
trend_basevalue <- reduce(trend_basevalue, `+`)

# standardize
trend_basevalue <- trend_basevalue / mean(trend_basevalue)

# ESTIMATE MU_DAILY
temp <- hist.weighted(day_marks, wi_daily, breaks = daily_base)
daily_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_daily)

daily_basevalue <- daily_basevalue / mean(daily_basevalue)


# ESTIMATE MU_WEEKLY
# as before, get intensity at points throughout the week
temp <- hist.weighted(week_marks, week_weights * wi_weekly, breaks = weekly_base)

# smooth density
weekly_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_weekly)
# standardize
weekly_basevalue <- weekly_basevalue / mean(weekly_basevalue)

# estimate mu_bg
background_basevalue <- matrix(0, 
                               nrow=length(background_base$x), 
                               ncol=length(background_base$y))

for(i in 1:nrow(a)){
  if (i %% 100 == 0) print(paste("on:", i))
  fn <- paste("Background.Smoothers/", "bgsmoother-", i, ".val", sep="")
  load(fn)
  background_basevalue <- background_basevalue + bg_probs[i] * bgsmoother
}

# Standardize the background so its average is 1 inside study area
background_basevalue <- background_basevalue/mean(background_basevalue[background_marks>=0])


# then update A and mu_0

# then calculate lambda_i = mu0 * background_i + A * trigger_i
# and lambda = mu0 * background_all + A * trigger_all

# terminate or loop back

# calculate g(t) and h(s) edge corrections
g_edge_correction <- unlist(map(a$days, function(x) sum(g_fun(seq(0, TT - x, time_units) + 0.6e-5)) * time_units))

h_edge_correction <- rep(0, nrow(a))

for(i in 1:nrow(a)){
  if(!file.exists("Excite.Spatail.Marks")){system("mkdir Excite.Spatail.Marks")}
  
  fn <- paste0("Excite.Spatail.Marks/crime1-", substr(10000 + i, 2, 5), ".mark")
    
  if(!file.exists(fn)){
    h_mark_temp <- matrix((inpoly(h_base_x + a$coorx[i], 
                                  h_base_y + a$coory[i], 
                                  city.boundary$x, 
                                  city.boundary$y) >=0),
                          ncol=length(h_base_x))
    save(h_mark_temp, file=fn)
  } else {
    load(fn)
  }
  h_edge_correction[i] <- simpson.2D(h_mark_temp * h_basevalue, space_units, space_units)
}

# get g(t) and h(s) propto equations

# evaluate triggering function
wi_gh <- A * g_fun(a$days[ij.mat[, 1]] - a$days[ij.mat[, 2]]) *
         h_fun(a$coorx[ij.mat[, 1]] - a$coorx[ij.mat[, 2]],
               a$coory[ij.mat[, 1]] - a$coory[ij.mat[, 2]]) /
         lambda_at_events[ij.mat[, 1]]

# I THINK THIS SHOULD BE G_EDGE_CORRECTION
temp <- hist.weighted(a$days[ij.mat[,1]] - a$days[ij.mat[,2]],
                      wi_gh / 
                      (h_edge_correction[ij.mat[, 2]] *
                       g_rep_fun(a$days[ij.mat[, 1]] - a$days[ij.mat[, 2]])),
                      breaks = g_base)

# where is this bw coming from?
g_basevalue <- ker.smooth.conv(temp$mids, 
                               temp$density, 
                               bandwidth=0.05)

g_basevalue <- g_basevalue/simpson(g_basevalue, time_units) 

# update of h_basevalue?
dis.mat <- cbind(a$coorx[ij.mat[,1]] - a$coorx[ij.mat[,2]], 
                 a$coory[ij.mat[,1]] - a$coory[ij.mat[,2]])

# I THINK THIS SHOULD BE SPATIAL.EDGE.CORRECTION
temp <- hist.weighted.2D(dis.mat[,1], 
                         dis.mat[,2], 
                         excite.wghs/(excite.temporal.edge.correction[ij.mat[,2]]), 
                         x.breaks= h_base_x, 
                         y.breaks= h_base_y)

# what is 2.35? something about the maximum distance
excite.spatial.mark2 <- (h_base_x^2 + h_base_y^2 < 2.35^2)

# expand repetition fun
temp <- h_rep_fun(temp$x.mids %o% rep(1, length(temp$y.mids)),
                  rep(1, length(temp$y.mids)) %o% temp$x.mids)

# cut off Kernel density estimates at distances of more than 2.35^2
h_basevalue <- ker.smooth.2D.fft(temp$x.mids,
                                 temp$y.mids, 
                                 temp$density, 
                                 x.bandwidth=0.1, 
                                 y.bandwidth=0.1) * excite.spatial.mark2

h_basevalue <- h_basevalue/simpson.2D(h_basevalue, space_units, space_units)


# get background_i and int_background
# Define a couple of functions which (linearly) interpolate between x and y
trend_fun <- approxfun(time_marks, trend_basevalue, yleft=0, yright=0)

weekly_fun <- function(x){
  approxfun(weekly_base, weekly_basevalue,             
            yleft=0, yright=0)(x - as.integer(x/7)*7)
}

daily_fun <- function(x){
  approxfun(daily_base, daily_basevalue,
            yleft=0, yright=0)(x - as.integer(x))
}

# interpolate mu_b
background_fun <- function(x,y) (interp.surface(obj=list(x = background_base$x, 
                                                         y = background_base$y,
                                                         z = background_basevalue),
                                                loc=cbind(x=c(x), y=c(y))))

# evaluate background at all events
bg_at_events_no_mu <- (trend_fun(a$days) * 
                         weekly_fun(a$days) * 
                         daily_fun(a$days) *
                         background_fun(a$coorx, a$coory))

# mean of background over entire space
bg_at_all_no_mu <- (mean(trend_fun(time_marks) * 
                           weekly_fun(time_marks) * 
                           daily_fun(time_marks)) * 
                      TT *
                      mean(background_fun(background_basex,
                                          background_basey) * 
                             background_marks) *
                      ra)

# get trigger_i and int_trigger

# interpolate g(t)
g_fun <- approxfun(seq(0, mac_t, time_units) + 0.6e-12,
                   g_basevalue, 
                   yleft=0, yright=0)

# interpolate h(s)
h_fun <- function(x,y){
  temp <- interp.surface(obj=list(x = h_base_x, 
                                  y = h_base_y, 
                                  z = h_basevalue),
                         loc=cbind(x=c(x), y=c(y))) 
  temp[is.na(temp)] <- 0
  temp
}

# trigger prob rho
rho_at_events_no_A <- map(1:nrow(a), 
                          function(i) g_fun(a$days - a$days[i]) * 
                            h_fun(a$coorx - a$coorx[i], 
                                  a$coory - a$coory[i]))
rho_at_events_no_A <- reduce(rho_at_events_no_A, `+`)

# this loop takes approx 3h because we need to evaluate 
# h_fun repeatedly over the entire grid
bg_weight <- sum(background_marks > 0)/length(background_marks)

# get mean effect of g_fun over the entire time
# multiply through with all constants
rho_at_all_no_A <- map(a$days, function(x) mean(g_fun(time_marks - x)) * TT * ra)

for (i in 1:nrow(a)){
  if (i %% 100 == 0) print(paste("on:", i))
  temp <- h_fun(background_basex - a$coorx[i], 
                background_basey - a$coory[i])
  rho_at_all_no_A[[i]] <- rho_at_all_no_A[[i]] * mean(temp[background_marks > 0]) * bg_weight
}

rho_at_all_no_A <- reduce(rho_at_all_no_A, sum)

# update mu0 and A according to equations (34) and (35)
old_mu0 <- mu0
old_A <- A

res <- optim(par=sqrt(c(mu0, A)), neg_loglik, control=list(trace=6))

# update
mu0 <- res$par[1]^2
A <- res$par[2]^2

lambda_at_events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
lambda_at_all <- mu0 * bg_at_all_no_mu + A * rho_at_all_no_A

bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events


if (abs(mu0 - old_mu0) < tol | abs(A - old_A) < tol){
  break
}

