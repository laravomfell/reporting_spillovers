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
voronoi_num_samples <- 5

initial_events_sf <- st_as_sf(da[da$e_type == 0,], coords = c("coorx", "coory"))

voronoi_polygons <- st_as_sf(st_intersection(st_cast(st_voronoi(st_union(initial_events_sf))), st_union(shp)))
voronoi_polygons$expected_count <- 0

# This will sample all cells in one go
#eval_points <- st_sample(voronoi_polygons, size=rep(voronoi_num_samples, nrow(voronoi_polygons)))

triggering_at_loc_with_all_times_fun <- function(a, x, y, i) {
    output <- g_fun(time_marks - a$days[i]) * h_fun(x - a$coorx[i], y - a$coory[i])
    return(output)
}
#temp <- triggering_at_loc_with_all_times_fun(da, 15.0, 15.0, 3)
#temp

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

for (cell_i in 1:5) {
    points_in_cell <- st_sample(voronoi_polygons[cell_i, ], size=voronoi_num_samples)

    # Evaluate the lambda at those locations.
    # TODO:
    
    # Average over the values of lambda.
    # TODO:

    # Store the value in the dataframe for the corresponding Voronoi cell.
    voronoi_polygons[cell_i, 'expected_count'] <- 1.0
}
plot(voronoi_polygons)




