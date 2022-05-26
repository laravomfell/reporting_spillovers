library(dplyr)
library(lubridate)
library(tidyr)


if (!file.exists("results/spatial_points_sf.rds")){
  
  # grab boundary
  region_boundary <- st_union(shp)
  #500 x 500 m cells
  grid_extent <- st_make_grid(shp, cellsize = 500)
  region_gridded <- st_as_sf(st_intersection(region_boundary, grid_extent))
  region_gridded$cell_id <- 1:nrow(region_gridded)
  # calculate area accounting for km^2 <-> m^2
  region_gridded$area <- st_area(region_gridded) 
  
  ## Turn the spatial intensity points into an SF object of points
  spatial_points_df <- expand.grid(x = background_base$x * 1000, 
                                   y = background_base$y * 1000)
  spatial_points_df <- spatial_points_df[as.vector(background_marks > 0), ]
  
  
  spatial_points_sf <- st_as_sf(spatial_points_df, coords = c("x", "y"), 
                                crs = st_crs(shp), agr = "constant") 
  
  spatial_points_cell_labels_sf <-st_join(spatial_points_sf, region_gridded,
                                          join = st_intersects, largest=T, left = T)
  
  st_write(region_gridded, paste0("results/exp_counts_spatial_grid.shp"),
           append = FALSE)
  saveRDS(spatial_points_cell_labels_sf, file = "results/spatial_points_sf.rds")
} else {
  spatial_points_cell_labels_sf <- readRDS("results/spatial_points_sf.rds")
}

# create temporal discretisation: chop off incomplete weeks from either side of the year.
x <- seq(as.Date(format.Date(da[1,]$datetime_unif, "%Y-01-01")), 
         as.Date(format.Date(da[1,]$datetime_unif, "%Y-12-31")), by = "day")
start_day <- which((weekdays(x) == "Monday" & yday(x) < 7) == TRUE)
end_day <- which((weekdays(x) == "Sunday" & yday(x) > length(x) - 7) == TRUE)
time_marks_cut <- seq(start_day - 1, end_day, 1/(TT/3.65))
coarse_times_weekly <- seq(start_day - 1, end_day, 7)
week_agg_labels <- cut(time_marks_cut, breaks = coarse_times_weekly, include.lowest = T, labels = F)

## Compute background components at all time_marks and all spatial points

# bg at all time marks
bg_temporal <- trend_fun(time_marks_cut) * weekly_fun(time_marks_cut) * daily_fun(time_marks_cut)

# bg at all grid points
bg_spatial_points <- background_fun(background_basex, background_basey)[as.vector(background_marks > 0)]

## Compute triggering components at all time_marks and all spatial points
## WARNING THIS IS MEMORY INTENSIVE: this will produce TT x SS matrix
# trigger_all_times_and_locs <- matrix(0, length(time_marks_cut), length(bg_spatial_points))
# for (i in 1:nrow(da)) {
#   if (i %% 500 == 0) print(paste("on:", i))
#   trigger_all_times_and_locs <- trigger_all_times_and_locs + theta[as.numeric(da[i,]$e_type) + 1] * trigger_general(a = da, i = i, time_points = time_marks_cut)
# }

## Combine to get lambda
lambda_all_times_and_locs <- mu0 * (bg_temporal %o% bg_spatial_points)# + trigger_all_times_and_locs

num_NAs <- sum(is.na(spatial_points_cell_labels_sf$cell_id))
print(paste("Number of unmatched points during the join: ", num_NAs))

## Project onto the coarser space-time grid with spatial join and temporal 
lambda_all_times_and_locs <- as.data.frame(t(lambda_all_times_and_locs))

lambda_all_times_and_locs$cell_id <- spatial_points_cell_labels_sf$cell_id

lambda_all_cells_and_times <- lambda_all_times_and_locs %>% drop_na() %>% group_by(cell_id) %>%
  summarise(across(everything(), mean))

# drop the 'cell_id' columns
lambda_all_cells_and_times <- subset(lambda_all_cells_and_times, select = -cell_id)

# multiply all rows by the corresponding cell area, dividing by 1000^2 (since our areas)
lambda_all_cells_and_times_area_int <- as.vector(spatial_points_cell_labels_sf$area) * lambda_all_cells_and_times / 1e06

# Transpose back to TT x SS shape
lambda_all_times_and_cells <- as.data.frame(t(lambda_all_cells_and_times_area_int))

# Aggregate into weeks
lambda_all_times_and_cells$week_id <- week_agg_labels  # 35304
lambda_all_weeks_and_cells <- lambda_all_times_and_cells %>% group_by(week_id) %>%
  summarise(across(everything(), mean))
lambda_all_weeks_and_cells <- 7 * lambda_all_weeks_and_cells
lambda_all_weeks_and_cells <- lambda_all_weeks_and_cells[ , !(names(lambda_all_weeks_and_cells) %in% c("week_id"))]

print(paste("Number of predicted events using the coarse discretisation buckets: ", sum(lambda_all_weeks_and_cells)))
print(paste("Number of actuall events in the dataset: ", sum(da$e_type == 0)))

# write as csv instead then
fwrite(lambda_all_weeks_and_cells, file = paste0("results/exp_counts_data_", experiment_id, ".csv"))
