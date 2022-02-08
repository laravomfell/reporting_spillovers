library(sf)
library(dplyr)
library(lubridate)


# shp <- read_sf("cov.shp", crs = 27700)



shp <- st_set_crs(shp, 27700)
region_boundary <- st_union(shp)
grid_extent <- st_make_grid(shp, cellsize = 0.5)
region_gridded <- st_as_sf(st_intersection(region_boundary, grid_extent))
region_gridded['cell_id'] <- 1:nrow(region_gridded)
region_gridded$area <- st_area(region_gridded)

# plot(region_gridded)


# create temporal discretisation: chop off incomplete weeks from either side of the year.
x <- seq(as.Date(format(da[1, 'date_unif'], "%Y-01-01")), as.Date(format(da[1, 'date_unif'], "%Y-12-31")), by = "day")
start_day <- which((weekdays(x) == "Monday" & yday(x) < 7) == TRUE)
end_day <- which((weekdays(x) == "Sunday" & yday(x) > length(x) - 7) == TRUE)
time_marks_cut <- seq(start_day - 1, end_day, 1/(TT/3.65))
coarse_times_weekly <- seq(start_day - 1, end_day, 7)
week_agg_labels <- cut(time_marks_cut, breaks = coarse_times_weekly, include.lowest = T, labels = F)



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


## Combine to get lambda
lambda_all_times_and_locs <- mu0 * (bg_temporal %o% bg_spatial_points) + trigger_all_times_and_locs




## Turn the spatial intensity points into an SF object of points
spatial_points_df <- expand.grid(x = background_base$x, 
                                 y = background_base$y)
spatial_points_df <- spatial_points_df[as.vector(background_marks > 0), ]


spatial_points_sf <- st_as_sf(spatial_points_df, coords = c("x", "y"), 
                                                crs = 27700, agr = "constant") 

spatial_points_cell_labels_sf <-st_join(spatial_points_sf, region_gridded,
                                        join = st_intersects) #%>% group_by(region_gridded.x) %>% summarise(mean(bg_spatial_points_sf.value))
spatial_points_cell_labels <- as.data.frame(spatial_points_cell_labels_sf$cell_id)
spatial_points_cell_areas <- as.data.frame(spatial_points_cell_labels_sf$area)


## Project onto the coarser space-time grid with spatial join and temporal 

lambda_all_locs_and_times <- transpose(lambda_all_times_and_locs)
lambda_all_locs_and_times['cell_id'] <- spatial_points_cell_labels
lambda_all_locs_and_times['area'] <- spatial_points_cell_areas

lambda_all_cells_and_times <- lambda_all_locs_and_times %>% group_by(cell_id) %>%
  summarise(across(everything(), mean))

# multiply all rows by the corresponding cell area
lambda_all_cells_and_times_area_int <- lambda_all_locs_and_times[, 'area'] * lambda_all_cells_and_times 

# drop the 'cell_id' and 'area' columns
lambda_all_cells_and_times_area_int <- lambda_all_cells_and_times_area_int[ , !(names(lambda_all_cells_and_times_area_int) %in% c("cell_id","area"))]

## Transpose back to TT x SS shape
lambda_all_times_and_cells <- transpose(lambda_all_cells_and_times_area_int)

# Aggregate into weeks
lambda_all_times_and_cells['week_id'] <- week_agg_labels
lambda_all_weeks_and_cells <- lambda_all_times_and_cells %>% group_by(week_id) %>%
  summarise(across(everything(), mean))
lambda_all_weeks_and_cells <- 7 *  lambda_all_weeks_and_cells
