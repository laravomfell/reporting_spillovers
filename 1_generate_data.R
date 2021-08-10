# In this file, we generate some fake data on a circle bounded to show the data 
# structure and how our model works


# helper function
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}


set.seed(42)

# in our data, each event has a unique case_id, but some events are follow-ups so they belong to another case
# in our data, the split is roughly 70/30
n <- num_all_events
n_initial <- floor(n * 0.7)

if (parents_proportion > 0.999) {
    print("Using homogeneous Poisson to generate the initial events.")
    cluster_process_sim <- rpp(1, s.region=as.matrix(boundary), t.region=c(0, end_date - start_date + 1),
                               npoints=n_initial, replace=TRUE, discrete.time=FALSE)
} else {
    print("Using a clustering process to generate the initial events.")
    lbda <- function(x,y,t,a){a*(2 + sin(((2 * pi) / 7) * t) + 2 + sin(2 * pi * t))}
    cluster_process_sim <- rpcp(s.region=as.matrix(boundary), t.region=c(0, end_date - start_date + 1),
                                nparents=floor(parents_proportion*n_initial), npoints=n_initial, infectious=TRUE,
                                cluster=c("normal", "exponential"),
                                dispersion=c(0.3, 10), lambda=lbda,
                                a=100000)  
}
# animation(cluster_process_sim$xyt, s.region=cluster_process_sim$s.region, t.region=cluster_process_sim$t.region, runtime=20)
# plot(cluster_process_sim$xyt, style="elegant")

initial <- data.table(
    grouped_id = 1:n_initial,
    e_type = 0L,
    coorx = cluster_process_sim$xyt[, 1],
    coory = cluster_process_sim$xyt[, 2],
    days =  cluster_process_sim$xyt[, 3]
)
initial[, datetime_unif := start_date + 3600 * 24 * days]



# create followups: sample events up to 14 dayes before the end of the period
cutoff_idx <- nrow(initial[datetime_unif < end_date - 3600 * 24 * 14])
followup <- initial[sample(cutoff_idx, size = floor(n - n_initial), replace = T),]
# add tiny amount of jitter to avoid duplicate locations
followup[, ":=" (coorx = coorx + runif(nrow(followup), -0.0005, 0.0005),
                 coory = coory + runif(nrow(followup), -0.0005, 0.0005))]
# add up between 1h to 14 days of time difference
followup[, datetime_unif := datetime_unif + runif(nrow(followup), 3600, 1209600)]

# change e_type to followup
followup[, e_type := 1L]


# Allow follow-up to create new first-time events (e_type=0).
triggered_by_follow_up <- followup[sample(nrow(followup),
                                          size = floor(followup_trig_prob * nrow(followup)),
                                          replace=T),]

# add Gaussian noise 
triggered_by_follow_up[, ":=" (coorx = coorx + rnorm(nrow(triggered_by_follow_up), 0, 0.15),
                               coory = coory + rnorm(nrow(triggered_by_follow_up), 0, 0.15))]

# add the exponentially distributed time increment with mean 10 days
triggered_by_follow_up[, datetime_unif := datetime_unif + rexp(nrow(triggered_by_follow_up), 10*24*3600)]

# cut out things outside the region and outside the time window
triggered_by_follow_up <- triggered_by_follow_up[datetime_unif <= end_date + 3600 * 24,]
points_sf <- st_as_sf(triggered_by_follow_up[, c('coorx', 'coory')], coords=c("coorx","coory"))
triggered_by_follow_up <- triggered_by_follow_up[as.numeric(st_intersects(points_sf, shp, sparse = FALSE)),]

# change e_type to triggered_by_follow_up
triggered_by_follow_up[, e_type := 0L]

da <- rbind(initial, followup, triggered_by_follow_up)
setorder(da, datetime_unif)
# create unique case id
da[, case_id := 1:.N]

# calculate event-specific bandwidth
# set lower bound on bandwidth
da$bandwidth <- rep(0.01, nrow(da))

print(paste("Using", n_p, "neighbours to define the smoothing disc."))
for(i in which(da$e_type == 0)){
  # calculate sq distances to all points which are initial events (including i)
  temp <- dist.squared(da$coorx[i], da$coory[i], 
                       da$coorx[da$e_type == 0], da$coory[da$e_type == 0])
  temp2 <- sqrt(sort(temp[temp >= 0.002^2])[n_p + 1]) # n_p +1'd  since we don't exclude i above
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
da$days <- (unclass(da$datetime_unif) - unclass(as.POSIXct(start_date)))/24/60/60

fwrite(da, file = paste0("da_", experiment_id, ".csv"), dateTimeAs = "write.csv")
