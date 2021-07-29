# This modified code adapted from Zhuang/Mateu proceeds as follows:

# INITALIZATION
# initialize mu_daily, mu_weekly, mu_trend, mu_bg
# initalize g(t) and h(s)

# get background: at events and integral
# get trigger: at events and integral
# then update mu_0 and theta based on loglik
# get lambda: at events and integral

# PRECALCULATE
# to later update g(t) and h(s) we need two quantities which we can precalculate now.
# 1. a precursor of the rho_ij matrix giving us "candidate" pairs of (i,j)
# which are close enough in space/time to be potentially triggered 
# 2. repetition terms calculating how often the trigger function is potentially "triggered"/ is applicable
# calculate \sum_i I(t_i + t \leq T) and \sum_i I(s_i + s \in s)


# ENTER LOOP:

# get weights wi_daily, wi_weekly, wi_trend, wi_bg
# use those to get mu_daily, mu_weekly, mu_trend, mu_bg
# get g(t) propto equation: calculate g(t) edge correction and get g(t)
# same for h(s)
# get background: at events and integral
# get trigger: at events and integral
# update mu_0 and theta

# terminate or loop back


# SETUP --------------------------------------------------------

event_types <- purrr::set_names(unique(da$e_type))
# since not all events are outcome events, we need a lookup vector
id <- 1:nrow(da)
id[da$e_type == 1] <- NA_integer_
id[!is.na(id)] <- 1:sum(!is.na(id))

# inpoly needs list to boundary coordinates
shp <- st_as_sf(w)
shp_area <- st_area(shp)

# extract boundaries
bbox <- st_bbox(shp)

boundary <- data.frame(st_coordinates(shp)[, c("X", "Y")])

# range of coordinates
x_range <- c(bbox["xmin"], bbox["xmax"])
y_range <- c(bbox["ymin"], bbox["ymax"])

# range of area
ra <- (x_range[2]-x_range[1])*(y_range[2]-y_range[1])

# TT = length of the time interval
TT <- ceiling(max(da$days))

# time stamps spaced out as 1/100 days (= 14min)
time_marks <- seq(0, TT, 1/(TT/3.65))

# maximum trigger ranges allowed
# 30 days and 2km
max_t <- 40
max_s <- 3

# at which resolution are we looking into triggering?
# 50m and ~7.5 min
time_units <- 0.005
space_units <- 0.05

# Initialization --------------------------------------------------------------

# 1. daily kernel

# daily_base is splitting the day into roughly 7.5min chunks
daily_base <- seq(0, 1, time_units)

# calculate how far into the day each event happened. unit: 1/days
day_marks <- da$days[da$e_type == 0] - as.integer(da$days[da$e_type == 0])

# get intensity at points throughout the day
temp <- hist.weighted(day_marks, rep(1, length(day_marks)), breaks = daily_base)

# smooth points from midpoints density with bw
daily_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_daily)

# standardize
daily_basevalue <- daily_basevalue/ mean(daily_basevalue)

# 2. weekly kernel

weekly_base <- seq(0, 7, time_units)

# calculate how far into the week each event happened. unit: 1/days
week_marks <- da$days[da$e_type == 0] - as.integer(da$days[da$e_type == 0] / 7) * 7 

week_weights <- 1 / (as.integer(TT/7) + 
                     (da$days[da$e_type == 0] - 
                     as.integer(da$days[da$e_type == 0] / 7) * 7 > TT - as.integer(TT / 7) * 7)) 

# as before, get intensity at points throughout the week
temp <- hist.weighted(week_marks, week_weights, breaks = weekly_base)

# smooth density
weekly_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_weekly)
# standardize
weekly_basevalue <- weekly_basevalue / mean(weekly_basevalue)

# 3. trend kernel

# calculate Gaussian kernel with boundary correction
# time marks is the variable 'x'. Here we are applying the kernel to prevent the leakage.
raw_trend_basevalue <- map(da$days[da$e_type == 0],
                           function(x) dnorm(x - time_marks, 0, bw_trend) + 
                                       dnorm(x + time_marks, 0, bw_trend) +
                                       dnorm(x + 2 * TT - time_marks, TT, bw_trend))

trend_basevalue <- reduce(raw_trend_basevalue, `+`)

# standardize
trend_basevalue <- trend_basevalue / mean(trend_basevalue)

# 4. area kernel

# define fine grid over entire study area
background_base <- list(x = seq(x_range[1], x_range[2], space_units), 
                        y = seq(y_range[1], y_range[2], space_units))

# create matrix of all locations on grid
background_basex <- background_base$x %o% rep(1, length(background_base$y))
background_basey <- rep(1, length(background_base$x))%o% background_base$y

# for each location, check if it's inside = 1, outside = 0 the boundary
background_marks <- Matrix(inpoly(background_basex,   
                                  background_basey, 
                                  boundary$X, 
                                  boundary$Y) >= 0, 
                           ncol = length(background_base$y),
                           sparse = TRUE)

# calculate background_basevalue
# (the initialization of the background kernel)

# predefine empty matrix (dgeMatrix is a standard dense matrix)
background_basevalue <- as(matrix(0,
                                  nrow = length(background_base$x), 
                                  ncol = length(background_base$y)),
                           "dgeMatrix")

# define parallelizable function
# write sparse matrix for each event
make_bg_smoothers <- function(i){
  # set file name
    fn <- paste0("background_smoothers/", gen_data_id, "_bgsmoother_", i, ".mtx")

    if(!file.exists(fn)) {
        ## for each location on the grid, add density of each event
        bgsmoother <- (dnorm(background_basex, 
                             da$coorx[i],
                             da$bandwidth[i]) * 
                       dnorm(background_basey,
                             da$coory[i], 
                             da$bandwidth[i]) /
                       da$bg_integral[i])
                                        # coerce values that are basically zero to be actually zero
        bgsmoother[bgsmoother < 0.001] <- 0
                                        # write as sparse matrix
        bgsmoother <- Matrix::Matrix(bgsmoother, sparse = T)
        Matrix::writeMM(bgsmoother, fn)
    }
}

# execute
foreach(i = which(da$e_type == 0)) %dopar% make_bg_smoothers(i)

# now read in each sparse matrix in a loop
for (i in which(da$e_type == 0)){
  if (i %% 500 == 0) print(paste("on:", i))
  
  fn <- paste0("background_smoothers/", gen_data_id, "_bgsmoother_", i, ".mtx")
  bgsmoother <- readMM(fn)
  
  background_basevalue <- background_basevalue + bgsmoother
}
  
# standardize the background so its average inside the study window is 1
background_basevalue <- background_basevalue / mean(background_basevalue[as.vector(background_marks > 0)])


# 5. initialize temporal triggering kernel
g_base <- seq(0, max_t, time_units)
g_basevalue <- (1/24 + g_base/24)^(-1)

# we normalize by the integral over the entire time window
g_basevalue <- g_basevalue / simpson(g_basevalue, time_units)

# 6. initalize spatial triggering
h_base_x <- seq(-max_s, max_s, space_units)
h_base_y <- seq(-max_s, max_s, space_units)
d <- length(h_base_x)

h_basevalue <- matrix(
  1/(abs(h_base_x %o% rep(1, d))^2 +
       abs(rep(1, d) %o% h_base_y)^2 + 1),
  ncol = d,
  nrow = d)

# normalize by integral over the entire spatial window
h_basevalue <- h_basevalue / simpson.2D(h_basevalue, space_units, space_units)


# load functions that will interpolate between these initialized values
source("interpolate_fun.R")

# 7. background

# evaluate background at all events (only the main events, not follow-ups)
bg_at_events_no_mu <- (trend_fun(da$days[da$e_type == 0]) * 
                       weekly_fun(da$days[da$e_type == 0]) * 
                       daily_fun(da$days[da$e_type == 0]) *
                       background_fun(da$coorx[da$e_type == 0], 
                                      da$coory[da$e_type == 0]))

# integral of background
bg_at_all_no_mu <- (mean(trend_fun(time_marks) * 
                         weekly_fun(time_marks) * 
                         daily_fun(time_marks)) * 
                    TT *
                    mean(background_fun(background_basex,
                                        background_basey) * 
                         background_marks) *
                    ra)

# 7. trigger

# at events, for both event types
rho_at_events_no_theta <- foreach(i = 1:nrow(da)) %dopar% trigger_fun(a = da, i = i)
  
# simplify
rho_at_events_no_theta <- map(event_types, 
                              function(x) reduce(rho_at_events_no_theta[da$e_type == x], `+`))
  
# integral of trigger (this operation takes approx 1.5h on our full dataset bc we need to evaluate
# h_fun repeatedly over the entire grid)

# precalculate some terms
bg_weight <- sum(background_marks > 0) / length(background_marks)

# get mean effect of g_fun and h_fun over the entire time and space
# multiply through with all constants
rho_at_all_no_theta <- foreach(i = 1:nrow(da)) %dopar% trigger_int_fun(da, time_marks, 
                                                                      background_basex, 
                                                                      background_basey, 
                                                                      background_marks, 
                                                                      TT * ra * bg_weight,
                                                                      i = i)
# reduce by event type
rho_at_all_no_theta <- unlist(map(event_types, function(x) reduce(rho_at_all_no_theta[da$e_type == x], sum)))

# Update mu0 and theta
mu0 <- 0.2
# just some initial values != 0
theta <- c(0.01, 0.01)

# the log likelihood function to minimize
neg_loglik <- function(x){
  # enforcing positivity of mu0 and theta
  mu0 <- x[1]^2
  theta <- x[-1]^2
  
  lambda_at_events <- mu0 * bg_at_events_no_mu + reduce(map2(rho_at_events_no_theta, theta, 
                                                             `*`),
                                                        `+`)
  
  lambda_at_all <- mu0 * bg_at_all_no_mu + reduce(map2(rho_at_all_no_theta, theta, 
                                                       `*`), 
                                                  `+`)
  
  - sum(log(lambda_at_events)) + lambda_at_all
}

# update mu0 and A
res <- optim(par = sqrt(c(mu0, theta)), neg_loglik)

mu0 <- res$par[1]^2
theta <- res$par[-1]^2


# get lambda
lambda_at_events <- mu0 * bg_at_events_no_mu + reduce(map2(rho_at_events_no_theta, theta, 
                                                           `*`), 
                                                      `+`)
lambda_at_all <- mu0 * bg_at_all_no_mu + reduce(map2(rho_at_all_no_theta, theta, 
                                                     `*`), 
                                                `+`)

# calculate bg_probs (phi_i)
bg_probs <- mu0 * bg_at_events_no_mu / lambda_at_events

# PRE-CALCULATION -------------------------------------------------------------

# 1. candidate (i,j) pairs

# create a data.frame ij where 
# `i` is the i index of all outcome events
# `j` lists all events j (including followups) where j < i
# and (i,j) are no more than 2km and 30 days apart
# this is a list of pairs where i could have been triggered by j
ij <- expand.grid(i = which(da$e_type == 0), j = 1:nrow(da))

ij <- ij[da$days[ij$i] > da$days[ij$j] &
         da$days[ij$i] < da$days[ij$j] + max_t &
         abs(da$coorx[ij$i] - da$coorx[ij$j]) < max_s &
         abs(da$coory[ij$i] - da$coory[ij$j]) < max_s,]

# get distance between (i,j)
dis <- data.frame(x = da$coorx[ij$i] - da$coorx[ij$j], 
                  y = da$coory[ij$i] - da$coory[ij$j])

# 2. repetition corrections

# repetance = "repetition corrections, i.e. for how many times the triggering effect at time lag t or the
# spatial offset (x, y) is observed"

# for each event time, count how many other event times are in the vicinity
g_rep <- reduce(map(da$days, function(x) as.numeric(g_base < TT - x)), `+`)
g_rep_fun <- approxfun(g_base, g_rep, yleft=1, yright=1)

# set up an empty matrix
h_rep <- matrix(0L, 
                ncol = length(h_base_x),
                nrow = length(h_base_y))
  
# this gives for each grid cell the number of events that can 'reach' this cell in terms of triggering.
for(i in 1:nrow(da)){
  fn <- paste0("h_space_marks/", experiment_id, "_h_marks_", i, ".csv")
  
  h_mark_temp <- matrix(inpoly(h_base_x %o% rep(1, d) + da$coorx[i],
                                 rep(1, d) %o% h_base_y + da$coory[i],
                                 boundary$X,
                                 boundary$Y) >= 0,
                          ncol = d)
  idx <- which(h_mark_temp == FALSE)
  # no need to write anything if all values are TRUE
  if (length(idx) == 0) next
  fwrite(data.table(idx), file = fn)

  h_rep <- h_rep + h_mark_temp
}


h_rep_fun <- function(x,y){
  temp <- interp.surface(obj=list(x = h_base_x, 
                                  y = h_base_y, 
                                  z = h_rep), 
                         loc=cbind(x=c(x), y=c(y)))
  temp[is.na(temp)] <- 0
  temp
}

h_rep_value <- h_rep_fun(da$coorx[ij$i] - da$coorx[ij$j], 
                         da$coory[ij$i] - da$coory[ij$j])


# a boolean map for the points whose distance is less than 2.35km from (0, 0).
excite.spatial.mark2 <- ((h_base_x %o% rep(1, d))^2 + (rep(1, d) %o% h_base_y)^2 < 2.35^2)

# INFERENCE LOOP --------------------------------------------------------------
k <- 1L

while (k < 40){
  print(paste(">>> Iteration:", k))
  print("updating background")
  
  # get weights wi_daily, wi_weekly, wi_trend, wi_bg
  # calculate w_i_t
  wi_trend <- trend_fun(da$days[da$e_type == 0]) * background_fun(da$coorx[da$e_type == 0], 
                                                                da$coory[da$e_type == 0]) / lambda_at_events

  # calculate w_i_d
  wi_daily <- daily_fun(da$days[da$e_type == 0]) * background_fun(da$coorx[da$e_type == 0], 
                                                                da$coory[da$e_type == 0]) / lambda_at_events
  
  # calculate w_i_w
  wi_weekly <- weekly_fun(da$days[da$e_type == 0]) * background_fun(da$coorx[da$e_type == 0], 
                                                                  da$coory[da$e_type == 0]) / lambda_at_events
  
  # use those to get mu_daily, mu_weekly, mu_trend, mu_bg
  
  # mu_daily: update intensity estimates with new weights
  temp <- hist.weighted(day_marks, wi_daily, breaks = daily_base)
  # smooth
  daily_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_daily)
  # standardize
  daily_basevalue <- daily_basevalue / mean(daily_basevalue)
  
  # mu_weekly: update intensity estimates with new weights
  temp <- hist.weighted(week_marks, week_weights * wi_weekly, breaks = weekly_base)
  # smooth density
  weekly_basevalue <- ker.smooth.fft(temp$mids, temp$density, bw_weekly)
  # standardize
  weekly_basevalue <- weekly_basevalue / mean(weekly_basevalue)

  # mu_trend: update intensity estimates with new weights
  trend_basevalue <- map2(raw_trend_basevalue, wi_trend, `*`)
  trend_basevalue <- reduce(trend_basevalue, `+`)
  # standardize
  trend_basevalue <- trend_basevalue / mean(trend_basevalue)
  
  # mu_area -- this recomputes the kernel after bg_probs has been updated
  background_basevalue <- as(matrix(0, 
                                    nrow = length(background_base$x), 
                                    ncol = length(background_base$y)),
                             "dgeMatrix")
  
  # read every single event density in again
  for(i in which(da$e_type == 0)){
    if (i %% 500 == 0) print(paste("on bgsmoother:", i))
    fn <- paste0("background_smoothers/", gen_data_id, "_bgsmoother_", i, ".mtx")
    bgsmoother <- readMM(fn)
    background_basevalue <- background_basevalue + bg_probs[id[i]] * bgsmoother
  }
  
  # standardize
  background_basevalue <- background_basevalue / mean(background_basevalue[as.vector(background_marks > 0)])
  
  
  print("calculating trigger edge corrections")
  
  # calculate g(t) and h(s) edge corrections
  g_edge_correction <- unlist(map(da$days, function(x) sum(g_fun(seq(0, TT - x, time_units) + 0.6e-5)) * time_units))
  
  
  calc_h_edge_corr <- function(i){
    # create matrix of 1's in correct dimensions
    fn <- paste0("h_space_marks/", experiment_id, "_h_marks_", i, ".csv")
    h_mark_temp <- matrix(1L, nrow = d, ncol = d)
    
    # if some elements are 0 (checked by file.exists),
    # set to 0
    if (file.exists(fn)){
      idx <- data.table::fread(fn)
      h_mark_temp[idx$idx] <- 0L
    }
    # and calculate Simpson integral approximation
    simpson.2D(h_mark_temp * h_basevalue, space_units, space_units)
  }
  
  # calculate edge correction in parallel
  h_edge_correction <- foreach(i = 1:nrow(da),
                               .combine = "c") %dopar% calc_h_edge_corr(i)
  
  print("updating trigger functions")
  # get g(t) and h(s) propto equations
  
  # evaluate triggering function * theta
  wi_gh <- theta[da$e_type[ij$j] + 1] * 
           g_fun(da$days[ij$i] - da$days[ij$j]) * 
           h_fun(da$coorx[ij$i] - da$coorx[ij$j],
                 da$coory[ij$i] - da$coory[ij$j]) /
           lambda_at_events[id[ij$i]]
  
  # correct wi_gh by g_edge and g_rep for g(t) update
  g_temp <- hist.weighted(da$days[ij$i] - da$days[ij$j],
                          wi_gh / 
                          (g_edge_correction[ij$j] *
                           g_rep_fun(da$days[ij$i] - da$days[ij$j])),
                          breaks = g_base)
  
  # update g(t)
  g_basevalue <- ker.smooth.conv(g_temp$mids, 
                                 g_temp$density, 
                                 bandwidth = bw_g)
  # standardize
  g_basevalue <- g_basevalue/simpson(g_basevalue, time_units) 
  
  
  # smooth wi_gh and correct for h_edge and h_rep
  # not originally in the code, but using edge + repetition correction
  temp <- hist.weighted.2D(dis$x, 
                           dis$y, 
                           wi_gh/(h_edge_correction[ij$j] * h_rep_value), 
                           x.breaks= h_base_x, 
                           y.breaks= h_base_y)
  
  # cut off Kernel density estimates at distances of more than 2.35^2
  h_basevalue <- ker.smooth.2D.fft(temp$x.mids,
                                   temp$y.mids, 
                                   temp$density, 
                                   x.bandwidth = bw_h, 
                                   y.bandwidth = bw_h) * excite.spatial.mark2
  # standardize
  h_basevalue <- h_basevalue/simpson.2D(h_basevalue, space_units, space_units)
  
  print("calculating background")  # with the updated background_basevalue, trend_basevalue, ...
  
  # get background_i and int_background
  # update the linear interpolation functions
  source("interpolate_fun.R")
  # evaluate background at all events
  bg_at_events_no_mu <- (trend_fun(da$days[da$e_type == 0]) * 
                         weekly_fun(da$days[da$e_type == 0]) * 
                         daily_fun(da$days[da$e_type == 0]) *
                         background_fun(da$coorx[da$e_type == 0], 
                                        da$coory[da$e_type == 0]))
  
  # integral of background over entire space
  bg_at_all_no_mu <- (mean(trend_fun(time_marks) * 
                             weekly_fun(time_marks) * 
                             daily_fun(time_marks)) * 
                        TT *
                        mean(background_fun(background_basex,
                                            background_basey) * 
                               background_marks) *
                        ra)
  
  print("calculating trigger") # withthe updated g_basevalue and h_basevalue
  # get trigger_i and int_trigger
  
  # trigger prob rho 
  # calculating using foreach
  rho_at_events_no_theta <- foreach(i = 1:nrow(da)) %dopar% trigger_fun(a = da, i = i)
  
  rho_at_events_no_theta <- map(event_types, function(x) reduce(rho_at_events_no_theta[da$e_type == x],
                                                                `+`))
  
  # get mean effect of g_fun and h_fun over the entire time and space
  # multiply through with all constants
  rho_at_all_no_theta <- foreach(i = 1:nrow(da)) %dopar% trigger_int_fun(da, time_marks, 
                                                                    background_basex, 
                                                                    background_basey, 
                                                                    background_marks, 
                                                                    TT * ra * bg_weight, 
                                                                    i = i)
  # reduce by event type
  rho_at_all_no_theta <- unlist(map(event_types, function(x) reduce(rho_at_all_no_theta[da$e_type == x], sum)))
  
  print("updating mu0 and theta")
  
  # update mu0 and theta
  old_mu0 <- mu0
  old_theta <- theta
  
  res <- optim(par = sqrt(c(mu0, theta)), neg_loglik)
  
  # update
  mu0 <- res$par[1]^2
  theta <- res$par[-1]^2
  
  print(paste("old mu0:", old_mu0, "new mu0:", mu0, "diff:", abs(mu0 - old_mu0)))
  print(paste("old theta0:", old_theta[1], "new theta0:", theta[1], "diff:", abs(theta[1] - old_theta[1])))
  print(paste("old theta1:", old_theta[2], "new theta0:", theta[2], "diff:", abs(theta[2] - old_theta[2])))
  
  # then calculate lambda_i = mu0 * background_i + A * trigger_i
  # and lambda = mu0 * background_all + A * trigger_all
  
  lambda_at_events <- mu0 * bg_at_events_no_mu + reduce(map2(rho_at_events_no_theta, theta, 
                                                             `*`), 
                                                        `+`)
  lambda_at_all <- mu0 * bg_at_all_no_mu + reduce(map2(rho_at_all_no_theta, theta, 
                                                       `*`), 
                                                  `+`)
  
  bg_probs <- mu0 * bg_at_events_no_mu / lambda_at_events
  
  # terminate or loop back
  conv <- abs(mu0 - old_mu0) < tol & all(abs(theta - old_theta) < tol)
  if (conv & k > 10){
    break
  }

  
  k <- k + 1L
  invisible(gc())
}
if (k == 41) print("loop ended after 40 iterations")

# keep track of influential events
influ_events <- foreach(i = 1:nrow(da)) %dopar% trigger_fun(a = da, i = i)
influ_events <- do.call(rbind, influ_events)


# Save results: mu0, theta, and the log-likelihood.
write(mu0, file = paste0("results/inferred_params", experiment_id, ".out"),
      append = FALSE)
write(theta, file = paste0("results/inferred_params", experiment_id, ".out"),
      append = TRUE, sep=" ", ncolumns=length(theta))
write(res$value, file = paste0("results/inferred_params", experiment_id, ".out"),
      append = TRUE)

if (save_snapshot) {
    save.image(file=paste0("results/snapshot_", experiment_id, ".RData"))
}
