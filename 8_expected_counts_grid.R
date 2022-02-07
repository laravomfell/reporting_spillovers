library(sf)
library(dplyr)
library(lubridate)


shp <- read_sf("cov.shp", crs = 27700)

region_boundary <- st_union(shp)
grid_extent <- st_make_grid(shp, cellsize = 500)
region_gridded <- st_as_sf(st_intersection(region_boundary, grid_extent))

plot(region_gridded)


# create temporal discretisation: chop off incomplete weeks from either side of the year.
x <- seq(as.Date(format(da[1, 'date_unif'], "%Y-01-01")), as.Date(format(da[1, 'date_unif'], "%Y-12-31")), by = "day")
start_day <- which((weekdays(x) == "Monday" & yday(x) < 7) == TRUE)
end_day <- which((weekdays(x) == "Sunday" & yday(x) > 358) == TRUE)
time_marks_cut <- seq(start_day, end_day, 1/(TT/3.65))


## Compute background components at all time_marks and all spatial points

# bg at all time marks
bg_temporal <-  trend_fun(time_marks_cut) * weekly_fun(time_marks) * daily_fun(time_marks)

# bg at all grid points
bg_spatial_points <-  background_fun(background_basex, background_basey)[as.vector(background_marks > 0)]


## Compute triggering components at all time_marks and all spatial points
## WARNING THIS IS REALLY SLOW
trigger_all_times_and_locs <- matrix(0, length(time_marks), length(bg_spatial_points))
for (i in 1:nrow(da)) {
  print(i)
  trigger_all_times_and_locs <- trigger_all_times_and_locs + theta[as.numeric(da[1, 'e_type'])] * trigger_general(a = da, i = i, time_points = time_marks)
}
# this will produce TT x SS matrix


lambda_all_times_and_locs <- (bg_temporal %o% bg_spatial_points) + trigger_all_times_and_locs



## Project onto the coarser space-time grid with spatial join and temporal 

# TODO: turn the spatial intensity points into an SF object of points
#intensity_spatial_points_sf <- st_as_sf(DT, coords = c("longitude", "latitude"), 
#                                                crs = 4326, agr = "constant") 
  
# TODO: spatial join into the grid, take the mean and multiplying by the cell area
# st_join(intensity_spatial_points_sf, region_gridded) %>% group_by(region_gridded.x) %>% summarise(mean(bg_spatial_points_sf.value))



## Combine to get lambda

# multiply through with theta and simplify
#lambda_at_all_locations <- mu0 * bg_at_all_locations + reduce(map2(trigger_at_all_locations, theta, `*`), `+`)

