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


# SETUP --------------------------------------------------------
library(spatstat) # must run R with option "--max-ppsize=100000"
library(purrr)
library(Matrix)
library(sf)
library(data.table)

library(foreach)
library(doParallel)
library(tictoc)

no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# parameter precision for mu0 and A
tol <- 0.001


# functions
source("utils.R")

tic("begin setup")
# read data (from source("prepare_data.R"))
a <- read.csv("da_cleaned.csv")
# Coventry boundary
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

# DEFINE PP PARS
# TT = length of the time interval
TT <- as.numeric(as.Date(max(a$datetime_unif)) - as.Date(min(a$datetime_unif))) + 1

# time stamps as 1 year/100 (= 3.65 days)
time_marks <- seq(0, TT, 1/(TT/3.65))

# range of coordinates
x_range <- c(bbox["xmin"], bbox["xmax"])
y_range <- c(bbox["ymin"], bbox["ymax"])

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
# changed from 0.002 originally
space_units <- 0.01


# STEP 1

# init mu_d
# 1. smoothing daily 

# daily_base is splitting the day into roughly 6min chunks
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
# We'll be reusing the raw, unnormalized quantities 
raw_trend_basevalue <- map(a$days, function(x) dnorm(x - time_marks, 0, 100) / (pnorm(TT, x, 100) - pnorm(0, x, 100)))
trend_basevalue <- reduce(raw_trend_basevalue, `+`)

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
# I've simplified this to (inside = 1, outside = 0) since no point lies on the
# boundary anyway
background_marks <- Matrix(inpoly(background_basex,   
                                  background_basey, 
                                  boundary$x, 
                                  boundary$y) >= 0, 
                           ncol = length(background_base$y),
                           sparse = TRUE)

# since the initialization of background_basevalue is slow, 
# check if I've already pre-calculated it

if (!file.exists("background_smoothers/setup_background_basevalue.csv")){
  # predefine empty matrix
  background_basevalue <- as(matrix(0,
                                    nrow = length(background_base$x), 
                                    ncol = length(background_base$y)),
                             "dgeMatrix")
  
  if(!file.exists("background_smoothers")){
    system("mkdir background_smoothers")
  }
  
  make_bg_smoothers <- function(i){
    fn <- paste0("background_smoothers/", "bgsmoother_", i, ".mtx")
    
    # for each location on the grid, add density of each event
    bgsmoother <- (dnorm(background_basex, 
                         a$coorx[i],
                         a$bandwidth[i]) * 
                     dnorm(background_basey,
                           a$coory[i], 
                           a$bandwidth[i]) /
                     a$bg_integral[i])
    # coerce values that are basically zero to be actually zero
    bgsmoother[bgsmoother < 0.001] <- 0
    bgsmoother <- as(bgsmoother, "sparseMatrix")
    writeMM(bgsmoother, fn)
  }
  
  smoothers <- foreach(i = 1:nrow(a)) %dopar% make_bg_smoothers(i)
  
  
  for(i in 1:nrow(a)){
    if (i %% 250 == 0) print(paste("on:", i))
    
    fn <- paste0("background_smoothers/", "bgsmoother_", i, ".mtx")
    bgsmoother <- readMM(fn)
    
    background_basevalue <- background_basevalue + bgsmoother
  }
  
  # standardize the background so its average inside the study window is 1
  background_basevalue <- background_basevalue / mean(background_basevalue[as.vector(background_marks > 0)])
  suppressMessages(fwrite(as.matrix(background_basevalue), 
                          file = "background_smoothers/setup_background_basevalue.csv"))
  
} else {
  background_basevalue <- as.matrix(fread("background_smoothers/setup_background_basevalue.csv"))
}



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
h_basevalue <- Matrix(
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
# calculating using foreach
if (!file.exists("setup_rho_at_events.Rdata")){
  rho_at_events_no_A <- foreach(i = 1:nrow(a)) %dopar% trigger_fun(a = a, i = i)
  rho_at_events_no_A <- reduce(rho_at_events_no_A, `+`)
  save(rho_at_events_no_A, file = "setup_rho_at_events.Rdata")
} else {
  load("setup_rho_at_events.Rdata")
}

# this loop takes approx 1.5h because we need to evaluate 
# h_fun repeatedly over the entire grid
bg_weight <- sum(background_marks > 0)/length(background_marks)


# get mean effect of g_fun and h_fun over the entire time and space
# multiply through with all constants
if (!file.exists("setup_rho_int.Rdata")){
  rho_at_all_no_A <- foreach(i = 1:nrow(a)) %dopar% trigger_int_fun(a, time_marks, 
                                                                    background_basex, 
                                                                    background_basey, 
                                                                    background_marks, 
                                                                    TT * ra * bg_weight, 
                                                                    i = i)
  rho_at_all_no_A <- reduce(rho_at_all_no_A, sum)
  save(rho_at_all_no_A, file = "setup_rho_int.Rdata")
} else {
  load("setup_rho_int.Rdata")
}


# Update A and mu for the first time
mu0 <- 0.77
A <- 0.03

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
# initial values really push down estimates of A
# maybe don't actually use loglik here and just take good initial values?
mu0 <- res$par[1]^2
A <- res$par[2]^2


# update other quantities
lambda_at_events <- mu0 * bg_at_events_no_mu + A * rho_at_events_no_A
lambda_at_all <- mu0 * bg_at_all_no_mu + A * rho_at_all_no_A

# calculate bg_probs (phi_i)
bg_probs <- mu0 * bg_at_events_no_mu / lambda_at_events


# pre-calculate the trigger terms

# candidate (i,j) pairs

# create a matrix ij where 
# `i` is the i index
# `j` lists all events j where j < i
# and (i,j) are no more than 2km and 15 days apart
temp <- (1:nrow(a)) %o% rep(1, nrow(a))

ij <- data.frame(i = c(t(temp)), j = c(temp))

ij <- ij[a$days[ij$i] > a$days[ij$j] &
           a$days[ij$i] < a$days[ij$j] + 15.0 &
           abs(a$coorx[ij$i] - a$coorx[ij$j]) < 2 &
           abs(a$coory[ij$i] - a$coory[ij$j]) < 2,]

# distance between (i,j)
dis <- data.frame(x = a$coorx[ij$i] - a$coorx[ij$j], 
                  y = a$coory[ij$i] - a$coory[ij$j])


# repetance = "repetition corrections, i.e. for how many times the triggering effect at time lag t or the
# spatial offset (x, y) is observed"

# for each event time, count how many other event times are in the vicinity
g_rep <- reduce(map(a$days, function(x) as.numeric(g_base < TT - x)), `+`)

# check if I've already precalculated h_rep previously
if (!file.exists("h_space_marks/setup_h_rep.csv")){
  
  h_rep <- matrix(0L, 
                  ncol = length(h_base_x),
                  nrow = length(h_base_y))
  
  if (!file.exists("h_space_marks")){
    system("mkdir h_space_marks")
  }
  
  for(i in 1:nrow(a)){
    if (i %% 250 == 0) print(paste("on:", i))
    
    fn <- paste0("h_space_marks/h_marks_", i, ".csv")
    
    h_mark_temp <- matrix(inpoly(h_base_x %o% rep(1, d) + a$coorx[i], 
                                 rep(1, d) %o% h_base_y + a$coory[i], 
                                 boundary$x, 
                                 boundary$y) >= 0,
                          ncol = d)
    idx <- which(h_mark_temp == FALSE)
    # no need to write anything if all values are TRUE
    if (length(idx) == 0) next
    fwrite(data.table(idx), file = fn)
    
    h_rep <- h_rep + h_mark_temp
  }
  
  suppressMessages(fwrite(h_rep, file = "h_space_marks/setup_h_rep.csv"))
} else {
  h_rep <- as.matrix(fread("h_space_marks/setup_h_rep.csv"))
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

h_rep_value <- h_rep_fun(a$coorx[ij$i] - a$coorx[ij$j], 
                         a$coory[ij$i] - a$coory[ij$j])

toc()

tic("enter the loop")
# ENTER THE LOOP ###############################################################
k <- 1L

while (k < 40){
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
  background_basevalue <- Matrix(0, 
                                 nrow=length(background_base$x), 
                                 ncol=length(background_base$y))
  
  for(i in 1:nrow(a)){
    if (i %% 100 == 0) print(paste("on:", i))
    fn <- paste0("background_smoothers/", "bgsmoother_", i, ".mtx")
    bgsmoother <- readMM(fn)
    background_basevalue <- background_basevalue + bg_probs[i] * bgsmoother
  }
  
  # Standardize the background so its average is 1 inside study area
  background_basevalue <- background_basevalue / mean(background_basevalue[background_marks > 0])
  
  
  # calculate g(t) and h(s) edge corrections
  g_edge_correction <- unlist(map(a$days, function(x) sum(g_fun(seq(0, TT - x, time_units) + 0.6e-5)) * time_units))
  
  h_edge_correction <- rep(0, nrow(a))
  
  for(i in 1:nrow(a)){
    fn <- paste0("h_space_marks/crime_", substr(10000 + i, 2, 5), ".mtx")
    h_mark_temp <- readMM(fn)
    h_edge_correction[i] <- simpson.2D(h_mark_temp * h_basevalue, space_units, space_units)
  }
  
  # get g(t) and h(s) propto equations
  
  # evaluate triggering function
  wi_gh <- A * g_fun(a$days[ij$i] - a$days[ij$j]) * 
    h_fun(a$coorx[ij$i] - a$coorx[ij$j],
          a$coory[ij$i] - a$coory[ij$j]) /
    lambda_at_events[ij$i]
  
  # I THINK THIS SHOULD BE G_EDGE_CORRECTION
  temp <- hist.weighted(a$days[ij$i] - a$days[ij$j],
                        wi_gh / 
                          (g_edge_correction[ij$j] *
                             g_rep_fun(a$days[ij$i] - a$days[ij$j])),
                        breaks = g_base)
  
  # where is this bw coming from?
  g_basevalue <- ker.smooth.conv(temp$mids, 
                                 temp$density, 
                                 bandwidth=0.05)
  
  g_basevalue <- g_basevalue/simpson(g_basevalue, time_units) 
  
  # update of h_basevalue
  # not originally in the code, but using edge + repetition correction
  temp <- hist.weighted.2D(dis$x, 
                           dis$y, 
                           wi_gh/(h_edge_correction[ij$j] * h_rep_value), 
                           x.breaks= h_base_x, 
                           y.breaks= h_base_y)
  
  # what is 2.35? something about the maximum distance
  excite.spatial.mark2 <- ((h_base_x %o% rep(1, d))^2 + (rep(1, d) %o% h_base_y)^2 < 2.35^2)

  # cut off Kernel density estimates at distances of more than 2.35^2
  h_basevalue <- ker.smooth.2D.fft(temp$x.mids,
                                   temp$y.mids, 
                                   temp$density, 
                                   x.bandwidth=0.1, 
                                   y.bandwidth=0.1) * excite.spatial.mark2
  
  h_basevalue <- h_basevalue/simpson.2D(h_basevalue, space_units, space_units)
  
  
  # get background_i and int_background
  # Define a couple of functions which (linearly) interpolate between x and y
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
  
  # get trigger_i and int_trigger
  
  # trigger prob rho 
  # calculating using foreach
  rho_at_events_no_A <- foreach(i = 1:nrow(a)) %dopar% trigger_fun(a = a, i = i)
  rho_at_events_no_A <- reduce(rho_at_events_no_A, `+`)
  
  # we need to evaluate 
  # h_fun repeatedly over the entire grid
  
  # get mean effect of g_fun and h_fun over the entire time and space
  # multiply through with all constants
  rho_at_all_no_A <- foreach(i = 1:nrow(a)) %dopar% trigger_int_fun(a, time_marks, 
                                                                    background_basex, 
                                                                    background_basey, 
                                                                    background_marks, 
                                                                    TT * ra * bg_weight, 
                                                                    i = i)
  rho_at_all_no_A <- reduce(rho_at_all_no_A, sum)
  
  
  # update mu0 and A according to equations (34) and (35)
  old_mu0 <- mu0
  old_A <- A
  
  res <- optim(par=sqrt(c(mu0, A)), neg_loglik, control=list(trace=6))
  
  # update
  mu0 <- res$par[1]^2
  A <- res$par[2]^2
  
  # then calculate lambda_i = mu0 * background_i + A * trigger_i
  # and lambda = mu0 * background_all + A * trigger_all
  
  lambda_at_events <- mu0 * bg_at_events_no_mu + A * rho_at_all_no_A
  lambda_at_all <- mu0 * bg_at_all_no_mu + A * rho_at_all_no_A
  
  bgprobs <- mu0 *  bg_at_events_no_mu / lambda_at_events
  
  # terminate or loop back
  if (abs(mu0 - old_mu0) < tol & abs(A - old_A) < tol){
    break
  }
  k <- k + 1L
}
toc()
if (k == 41) print("loop ended after 40 iterations")




