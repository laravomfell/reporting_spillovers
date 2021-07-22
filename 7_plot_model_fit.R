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

#save.image(file='inferred_model__stat.RData')
load('inferred_model__stat.RData')

voronoi_num_samples <- 5
initial_events_sf <- st_as_sf(da[da$e_type == 0,], coords = c("coorx", "coory"))
voronoi_polygons <- st_as_sf(st_intersection(st_cast(st_voronoi(st_union(initial_events_sf))), st_union(shp)))

voronoi_polygons$poly_id <- 1:nrow(voronoi_polygons)
voronoi_polygons$area <- st_area(voronoi_polygons)


#
# Slow version
#
voronoi_polygons$expected_count <- 0
triggering_at_loc_with_all_times_fun <- function(a, x, y, i) {
    output <- g_fun(time_marks - a$days[i]) * h_fun(x - a$coorx[i], y - a$coory[i])
    return(output)
}

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

for (cell_i in 1:nrow(voronoi_polygons) {
    points_in_cell <- st_sample(voronoi_polygons[cell_i, ], size=voronoi_num_samples)

    polygon_area <- st_area(voronoi_polygons[cell_i, ])

    # Evaluate the lambda at those locations.
    lambda_evals <- rep(0.0, voronoi_num_samples)
    for (i in 1:voronoi_num_samples) {
        point_coords <- st_coordinates(points_in_cell[i,])
        lambda_evals[i] <- lambda_at_loc(point_coords[1], point_coords[2])
    }    

    # Store the value in the dataframe for the corresponding Voronoi cell.
    voronoi_polygons[cell_i, 'expected_count'] <- polygon_area * mean(lambda_evals)
}
plot(voronoi_polygons)



#
# Faster version
#
triggering_at_locs_with_all_times_fun <- function(a, xs, ys, i) {
    h_all_locs <- h_fun(xs - a$coorx[i], ys - a$coory[i])
    g_at_time_marks <- g_fun(time_marks - a$days[i])
    trigger_constant <- theta[unlist(a[i, 'e_type']) + 1]
    output <- colSums(trigger_constant * (g_at_time_marks %o% h_all_locs)) # gives a T x S matrix
    return(output)
}

# This will sample all cells in one go
eval_points <- st_sample(voronoi_polygons, size=rep(voronoi_num_samples, nrow(voronoi_polygons)))
eval_points_coords <- st_coordinates(eval_points)

time_integral_bg <- mean(trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks)) * TT


pb <- txtProgressBar(min=1, max=nrow(da), initial=1)
trig_part_all_times <- Matrix(0, 1, nrow(eval_points_coords))
for (i in 1:nrow(da)) {
    trig_part_all_times <- trig_part_all_times + triggering_at_locs_with_all_times_fun(a = da, x = eval_points_coords[, 1], y = eval_points_coords[, 2], i = i)
    setTxtProgressBar(pb, i)
}

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

voronoi_expected_counts_sf <- merge(voronoi_polygons, polygon_expected_count)

voronoi_expected_counts_sf$residuals <- 1 - voronoi_expected_counts_sf$expected_count


pdf(paste0("figures/", experiment_id, "_voronoi_residuals.pdf")) 
plot(voronoi_expected_counts_sf["residuals"], )
dev.off()





theta = seq(0,3,length=500)


p <- ggplot(voronoi_expected_counts_sf, aes(x=residuals)) + 
    #geom_histogram(color="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") +
    stat_function(aes(y = ..y..), fun = dgamma, n = 101, args = list(shape = 3.659, scale = 1/3.659), alpha=0.4)
p
