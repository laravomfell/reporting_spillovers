# In this file, we perform the model check shown in figure 2.4

bg_at_all_locations <- trend_fun(time_marks) * 
  weekly_fun(time_marks) * 
  daily_fun(time_marks) * 
  mean(background_fun(background_basex, 
                      background_basey)[as.vector(background_marks > 0)]) * 
  ra * bg_weight


trigger_at_all_locations <- foreach(i = 1:nrow(da)) %dopar% trigger_at_all_fun(a = da, i = i, constants = bg_weight * ra)
# reduce by event type
trigger_at_all_locations <- map(event_types, function(x) reduce(trigger_at_all_locations[da$e_type == x],
                                                                `+`))

# multiply through with theta and simplify
lambda_at_all_locations <- mu0 * bg_at_all_locations + reduce(map2(trigger_at_all_locations, theta, `*`), `+`)

# creating the diagonal plot
n_events <- nrow(da[da$e_type == 0,])

fit <- data.frame(y = cumsum(lambda_at_all_locations) * (time_marks[2] - time_marks[1]),
                  x = stepfun(da$days[da$e_type == 0], 0:n_events)(time_marks))
p <- ggplot(fit,
            aes(x,y)) + 
  geom_ribbon(data = data.frame(x = 0:n_events, y = 0,
                                ymax = n_events * qbeta(.975, 0:n_events + 1, n_events - (0:n_events) + 1),
                                ymin = n_events * qbeta(.025, 0:n_events + 1, n_events - (0:n_events) + 1)),
              aes(x, ymin = ymin, ymax = ymax), fill = "#d6d6d6") +
  geom_line(data = data.frame(x = c(0, n_events),
                              y = c(0, n_events)),
            color = "black") +
  geom_line(color = "#721F81FF", size = .8) + 
  #scale_x_continuous(breaks = c(0, 2000, 4000, n_events)) + 
  #scale_y_continuous(breaks = c(0, 2000, 4000, n_events)) +
  coord_fixed(xlim = c(0, max(fit$x)), ylim = c(0, n_events))+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9)) +
  labs(x = "Event index", y = "Cumulative frequency")

ggsave(paste0("figures/", experiment_id, "_trans_time.pdf"), plot = p, width = 3.5, height = 3.5)

# plotting the deviation from the diagonal
conf <- data.frame(x = 0:n_events, y = 0,
                   ymax = n_events * qbeta(.975, 0:n_events + 1, n_events - (0:n_events) + 1),
                   ymin = n_events * qbeta(.025, 0:n_events + 1, n_events - (0:n_events) + 1))
fit$dev <- fit$y - fit$x
p <- ggplot(fit, aes(y, dev)) + 
  geom_ribbon(data = conf,
              aes(x, ymin = ymin -x, ymax = ymax -x, y=y),inherit.aes = F, fill = "#d6d6d6") +
  geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = c(0, 1000, 2000, n_events)) + 
  scale_y_continuous(breaks = c(-50,-25,0,25,50)) +
  geom_line(color = "#721F81FF") +
  labs(x = "Event index", y = parse(text = "Predicted~-~expected~number~of~events"))

ggsave(paste0("figures/", experiment_id, "_dev.pdf"), plot = p, width = 5.5, height = 3.5)


##
##
## Voronoi diagnostics
##
##
if (compute_voronoi) {
  voronoi_num_samples <- 50
  initial_events_sf <- st_as_sf(da[da$e_type == 0,], coords = c("coorx", "coory"))
  voronoi_polygons <- st_as_sf(st_intersection(st_cast(st_voronoi(st_union(initial_events_sf))), st_union(shp)))
  
  voronoi_polygons$poly_id <- 1:nrow(voronoi_polygons)
  voronoi_polygons$area <- st_area(voronoi_polygons)
  
  # 1- slow, 2- fast but memory-intensive, 3 (or anything else)- quite fast but memory-efficient
  voronoi_computation_variant <- 3
  
  
  time_integral_bg <- mean(trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks)) * TT
  
  
  # Computes the triggering (without theta) for all integration time marks
  triggering_at_loc_with_all_times_fun <- function(a, x, y, i) {
    output <- g_fun(time_marks - a$days[i]) * h_fun(x - a$coorx[i], y - a$coory[i])
    return(output)
  }
  
  # Computer triggering with corresponding theta coefficient for all
  # locations and averaged over the integration time marks
  triggering_at_locs_with_all_times_fun <- function(a, xs, ys, i, theta) {
    h_all_locs <- h_fun(xs - a$coorx[i], ys - a$coory[i])
    g_at_time_marks <- g_fun(time_marks - a$days[i])
    trigger_constant <- theta[unlist(a[i, 'e_type']) + 1]
    output <- colSums(trigger_constant * (g_at_time_marks %o% h_all_locs)) # gives a T x S matrix
    return(output)
  }
  
  # Computes the the intensity at a specific location (x, y), integrated
  # over the time domain.
  lambda_at_loc <- function(x, y) {
    time_integral_bg <- mean(trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks)) * TT
    
    trig_part_no_theta <- foreach(i = 1:nrow(da), .export=ls(envir=globalenv())) %dopar% {
      triggering_at_loc_with_all_times_fun(a = da, x = x, y = y, i = i)
    }
    
    trig_part_no_theta <- map(event_types,
                              function(x) reduce(trig_part_no_theta[da$e_type == x], `+`))
    output <- mu0 * background_fun(x, y) * time_integral_bg + TT * mean(reduce(map2(trig_part_no_theta, theta, `*`), `+`))
    return(output)
  }
  
  if (voronoi_computation_variant == 1) { # slow iterative version
    print("Computing Voronoi expected counts using the slow iterative method.")
    voronoi_polygons$expected_count <- 0
    for (cell_i in 1:nrow(voronoi_polygons)) {
      points_in_cell <- st_sample(voronoi_polygons[cell_i, ], size=voronoi_num_samples)
      
      polygon_area <- voronoi_polygons[cell_i, 'area']
      
      # Evaluate the lambda at those locations.
      lambda_evals <- rep(0.0, voronoi_num_samples)
      for (i in 1:voronoi_num_samples) {
        point_coords <- st_coordinates(points_in_cell[i,])
        lambda_evals[i] <- lambda_at_loc(point_coords[1], point_coords[2])
      }    
      
      # Store the value in the dataframe for the corresponding Voronoi cell.
      voronoi_polygons[cell_i, 'expected_count'] <- polygon_area * mean(lambda_evals)
    }
  } else if (voronoi_computation_variant == 2) {
    print("Computing Voronoi expected counts using the fast memory-intensive method.")
    
    # This will sample all cells in one go
    eval_points <- st_sample(voronoi_polygons, size=rep(voronoi_num_samples, nrow(voronoi_polygons)))
    eval_points_coords <- st_coordinates(eval_points)
    
    
    
    
    pb <- txtProgressBar(min=1, max=nrow(da), initial=1)
    trig_part_all_times <- Matrix(0, 1, nrow(eval_points_coords))
    for (i in 1:nrow(da)) {
      trig_part_all_times <- trig_part_all_times + triggering_at_locs_with_all_times_fun(a = da,
                                                                                         x = eval_points_coords[, 1],
                                                                                         y = eval_points_coords[, 2],
                                                                                         i = i, theta = theta)
      setTxtProgressBar(pb, i)
    }
    
    # If the computer with lots of memory is available then this
    # could be run in parallel, but note each of the subprocesses
    # creates a matrix of size T x S where S is the number of
    # *ALL* evaluation points which can be quite large as there is
    # `voronoi_num_samples` of evaluation for *each* cell.
    
    #trig_part_all_times <- foreach(i = 1:nrow(da),
    #                              .export=c('h_fun', 'g_fun', 'time_marks', 'theta'),
    #                              .packages=c('Matrix'),
    #                              .combine='+') %dopar% {
    #    triggering_at_locs_with_all_times_fun(a = da, x = eval_points_coords[, 1], y = eval_points_coords[, 2], i = i)
    #}
    
    trigger_contribution <- (TT / length(time_marks)) * trig_part_all_times # This should leave a vector of size S
    lambda_at_probes <- mu0 * background_fun(eval_points_coords[, 1], eval_points_coords[, 2]) * time_integral_bg + trigger_contribution
    
    intensity_sf <- st_sf(intensity=lambda_at_probes, eval_points)
    intensity_with_blocks_sf <- st_join(intensity_sf, voronoi_polygons, join = st_within)
    polygon_expected_count <- data.table(intensity_with_blocks_sf)[, list(expected_count=mean(intensity) * head(area, 1)), by='poly_id']
    
    voronoi_polygons <- merge(voronoi_polygons, polygon_expected_count)
  } else {
    print("Computing Voronoi expected counts using memory efficient method.")
    
    # This version will parallelise over the cells. It has low
    # memory footprint compared to the version above.  for each
    # cell we sample S points, create a T x S matrix and compute
    # the lambdas there.
    tic()
    voronoi_expected_counts <- foreach(cell_i= 1:nrow(voronoi_polygons),
                                       .packages=c('Matrix', 'sf', 'fields')) %dopar% {
       points_in_cell <- st_sample(voronoi_polygons[cell_i, ], size=voronoi_num_samples)
       eval_points_coords <- st_coordinates(points_in_cell)
       polygon_area <- st_drop_geometry(voronoi_polygons[cell_i, 'area'])
       triggering <- numeric(nrow(eval_points_coords))
       for (i in 1:nrow(da)) {
         triggering <- triggering + triggering_at_locs_with_all_times_fun(a = da,
                                                                          x = eval_points_coords[, 1], y = eval_points_coords[, 2],
                                                                          i = i, theta = theta)
       }
       trigger_contribution <- (TT / length(time_marks)) * triggering
       lambda_at_probes <- mu0 * background_fun(eval_points_coords[, 1],
                                                eval_points_coords[, 2]) * time_integral_bg + trigger_contribution
       expected_count_for_cell <-  polygon_area * mean(lambda_at_probes)
       expected_count_for_cell
    }
    toc()
    voronoi_polygons[, 'expected_count'] <- unlist(voronoi_expected_counts)
  }
  
  
  # Compute the residuals and make plots
  voronoi_polygons$residuals <- 1 - voronoi_polygons$expected_count
  
  pdf(paste0("figures/", experiment_id, "_voronoi_residuals.pdf")) 
  plot(voronoi_polygons["residuals"], )
  dev.off()
  
  pdf(paste0("figures/", experiment_id, "_voronoi_expected_gamma.pdf")) 
  p <- ggplot(voronoi_polygons) +
    geom_histogram(aes(x = residuals, y = ..density..),
                   binwidth = 0.1, fill = "grey", color = "black") +
    stat_function(aes(x = residuals, y = ..y..),
                  fun = function(x) dgamma(1-x, shape=3.569, scale=1/3.569),
                  n = 100, alpha=0.8, color = "red")
  p
  dev.off()
  
  st_write(voronoi_polygons, paste0("results/voronoi_", experiment_id, ".shp"))
}
