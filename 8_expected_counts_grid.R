library(sf)
library(dplyr)
library(lubridate)
library(tidyr)


shp <- st_set_crs(shp, 27700)
region_boundary <- st_union(shp)
grid_extent <- st_make_grid(shp, cellsize = 0.5)
region_gridded <- st_as_sf(st_intersection(region_boundary, grid_extent))
region_gridded$cell_id <- 1:nrow(region_gridded)
region_gridded$area <- st_area(region_gridded)


# create temporal discretisation: chop off incomplete weeks from either side of the year.
x <- seq(as.Date(format(da[1, 'date_unif'], "%Y-01-01")), as.Date(format(da[1, 'date_unif'], "%Y-12-31")), by = "day")
start_day <- which((weekdays(x) == "Monday" & yday(x) < 7) == TRUE)
end_day <- which((weekdays(x) == "Sunday" & yday(x) > length(x) - 7) == TRUE)
time_marks_cut <- seq(start_day - 1, end_day, 1/(TT/3.65))
coarse_times_weekly <- seq(start_day - 1, end_day, 7)
week_agg_labels <- cut(time_marks_cut, breaks = coarse_times_weekly, include.lowest = T, labels = F)



## Compute background components at all time_marks and all spatial points

# bg at all time marks
bg_temporal <-  trend_fun(time_marks_cut) * weekly_fun(time_marks_cut) * daily_fun(time_marks_cut)

# bg at all grid points
bg_spatial_points <-  background_fun(background_basex, background_basey)[as.vector(background_marks > 0)]


## Compute triggering components at all time_marks and all spatial points
## WARNING THIS IS MEMORY INTENSIVE: this will produce TT x SS matrix
trigger_all_times_and_locs <- matrix(0, length(time_marks_cut), length(bg_spatial_points))
for (i in 1:nrow(da)) {
  print(i)
  trigger_all_times_and_locs <- trigger_all_times_and_locs + theta[as.numeric(da[i, 'e_type'])] * trigger_general(a = da, i = i, time_points = time_marks_cut)
}

## Combine to get lambda
lambda_all_times_and_locs <- mu0 * (bg_temporal %o% bg_spatial_points) + trigger_all_times_and_locs


## Turn the spatial intensity points into an SF object of points
spatial_points_df <- expand.grid(x = background_base$x, 
                                 y = background_base$y)
spatial_points_df <- spatial_points_df[as.vector(background_marks > 0), ]


spatial_points_sf <- st_as_sf(spatial_points_df, coords = c("x", "y"), 
                              crs = 27700, agr = "constant") 

spatial_points_cell_labels_sf <-st_join(spatial_points_sf, region_gridded,
                                        join = st_intersects, largest=T, left = T)

num_NAs <- sum(is.na(spatial_points_cell_labels_sf$cell_id))
print(paste("Number of unmatched points during the join: ", num_NAs))

spatial_points_cell_labels <- spatial_points_cell_labels_sf$cell_id
spatial_points_cell_areas <- spatial_points_cell_labels_sf$area


## Project onto the coarser space-time grid with spatial join and temporal 
lambda_all_locs_and_times <- t(lambda_all_times_and_locs)

lambda_all_locs_and_times <- as.data.frame(lambda_all_locs_and_times)
lambda_all_locs_and_times$cell_id <- spatial_points_cell_labels
lambda_all_locs_and_times$area <- spatial_points_cell_areas

lambda_all_cells_and_times <- lambda_all_locs_and_times %>% drop_na() %>% group_by(cell_id) %>%
  summarise(across(everything(), mean))

# drop the 'cell_id' and 'area' columns after saving the areas
cell_areas <- as.numeric(lambda_all_cells_and_times$area)
lambda_all_cells_and_times <- lambda_all_cells_and_times[ , !(names(lambda_all_cells_and_times) %in% c("cell_id","area"))]

# multiply all rows by the corresponding cell area
lambda_all_cells_and_times_area_int <-  cell_areas * lambda_all_cells_and_times

# Transpose back to TT x SS shape
lambda_all_times_and_cells <- as.data.frame(t(lambda_all_cells_and_times_area_int))

# Aggregate into weeks
lambda_all_times_and_cells$week_id <-week_agg_labels  # 35304
lambda_all_weeks_and_cells <- lambda_all_times_and_cells %>% group_by(week_id) %>%
  summarise(across(everything(), mean))
lambda_all_weeks_and_cells <- 7 *  lambda_all_weeks_and_cells
lambda_all_weeks_and_cells <- lambda_all_weeks_and_cells[ , !(names(lambda_all_weeks_and_cells) %in% c("week_id"))]

print(paste("Number of predicted events using the coarse discretisation buckets: ", sum(lambda_all_weeks_and_cells)))
print(paste("Number of actuall events in the dataset: ", sum(da$e_type == 0)))

st_write(region_gridded, paste0("results/exp_counts_spatial_grid_", experiment_id, ".shp"))
save(lambda_all_weeks_and_cells, file = paste0("results/exp_counts_data", experiment_id, ".RData"))
