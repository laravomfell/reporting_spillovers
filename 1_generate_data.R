# In this file, we generate some fake data on a circle bounded by [0,1] to show the data 
# structure and how our model works

# function to generate random dates adapted from 
# https://stackoverflow.com/questions/14720983/efficiently-generate-a-random-sample-of-times-and-dates-between-two-dates
gen_date <- function(N, start = "2018-01-01", end = "2018-12-31"){
  st <- as.POSIXct(as.Date(start))
  et <- as.POSIXct(as.Date(end))
  # calculate distance between start and end
  dt <- as.numeric(difftime(et, st, unit = "sec"))
  # random draw over distance
  ev <- sort(runif(N, 0, dt))
  # return ordered dates
  st + ev
}

# this function already assumes that we're generating pts
# on a circle bounded by [0,1] with radius 0.5
gen_pts <- function(N){
  # pick radius of pts as function of circle radius 0.5
  rad <- 0.5 * sqrt(runif(N))
  # pick angles
  ang <- runif(N) * 2 * pi
  # create coords
  x <- 0.5 + rad * cos(ang)
  y <- 0.5 + rad * sin(ang)
  return(list(x = x, y = y))
}

# create circle with radius 0.5 on [0,1]
angle_increments <- 2 * pi/100
# create angles
ang <- seq(0, 2 * pi - angle_increments, by = angle_increments)
circ <- data.frame(x = 0.5 + 0.5 * cos(ang),
                   y = 0.5 + 0.5 * sin(ang))

# pretty self explanatory
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}


set.seed(1)

# in our data, each event has a unique case_id, but some events are follow-ups so they belong to another case
# in our data, the split is roughly 70/30
n <- 500
n_initial <- floor(n * 0.7)

# generate initial events
initial <- data.table(
  grouped_id = 1:n_initial,
  datetime_unif = gen_date(n_initial),
  e_type = 0L
)

initial[, c("coorx", "coory") := gen_pts(n_initial)]

# create followups
followup <- initial[sample(n_initial, size = floor(n - n_initial), replace = T),]
# add tiny amount of jitter to avoid duplicate locations
followup[, ":=" (coorx = coorx + runif(nrow(followup), -0.0005, 0.0005),
                 coory = coory + runif(nrow(followup), -0.0005, 0.0005))]
# add up between 1h to 14 days of time difference
followup[, datetime_unif := datetime_unif + runif(nrow(followup), 3600, 1209600)]
# but cap at 2018-12-31
followup[datetime_unif >= "2019-01-01", datetime_unif := as.POSIXct("2018-12-31 00:00:00") + runif(.N, 0, 86400)]
# change e_type to followup
followup[, e_type := 1L]

da <- rbind(initial, followup)
setorder(da, datetime_unif)
# create unique case id
da[, case_id := 1:.N]


# now prepare the data
# here we would normally load our shapefile, but instead we use the generated circle
w <- owin(poly = circ)

# calculate event-specific bandwidth
# set lower bound on bandwidth
da$bandwidth <- rep(0.01, nrow(da))

for(i in which(da$e_type == 0)){
  # calculate sq distances to all points which are initial events (including i)
  temp <- dist.squared(da$coorx[i], da$coory[i], 
                       da$coorx[da$e_type == 0], da$coory[da$e_type == 0])
  # n_p = 5 but need to use six since we don't exclude i above
  temp2 <- sqrt(sort(temp[temp >= 0.002^2])[6])
  # replace bandwidth if too small otherwise
  if(da$bandwidth[i] <= temp2){
    da$bandwidth[i] <- temp2
  }
}

# next, calculate integral around each event. this takes a while
da$bg_integral <- rep(NA_real_, nrow(da))

# calculate exact integral
for(i in which(da$e_type == 0)){
  # calculate exact integral
  da$bg_integral[i] <- polyCub.exact.Gauss(w,
                                           mean=c(da$coorx[i], da$coory[i]),
                                           Sigma=diag(da$bandwidth[i], 2),
                                           plot=F)
}

# transform datestamp into 1/days units since beginning of study window
da$days <- (unclass(da$datetime_unif) - unclass(min(as.POSIXct(as.Date(da$datetime_unif)))))/24/60/60

fwrite(da, file = "da.csv", dateTimeAs = "write.csv")

