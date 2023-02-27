## WARNING:
## This is very memory intensive and very much depends on the discretisation of 
## the spatial domain mostly. One can use the default background_basex/background_basey 
## if computer memory allows that. Otherwise, it's ok to use something less precise.
## For example to get discretisation half the resolution of the model's spatial
## resolution one can do as follows:

# define fine grid over entire study area
COARSE_CONSTANT <- 1.0
ripley_bg_base <- list(x = seq(x_range[1], x_range[2], COARSE_CONSTANT * space_units), 
                       y = seq(y_range[1], y_range[2], COARSE_CONSTANT * space_units))

# create matrix of all locations on grid
ripley_bg_basex <- ripley_bg_base$x %o% rep(1, length(ripley_bg_base$y))
ripley_bg_basey <- rep(1, length(ripley_bg_base$x)) %o% ripley_bg_base$y

# for each location, check if it's inside = 1, outside = 0 the boundary
ripley_bg_marks <- Matrix(inpoly(ripley_bg_basex,   
                                 ripley_bg_basey, 
                                 boundary$X, 
                                 boundary$Y) >= 0, 
                          ncol = length(ripley_bg_base$y),
                          sparse = TRUE)

# Compute triggering with corresponding theta coefficient for all
# locations and averaged over the integration time marks
triggering_at_locs_with_all_times_fun <- function(a, xs, ys, i, theta) {
  h_all_locs <- h_fun(xs - a$coorx[i], ys - a$coory[i])
  g_at_time_marks <- g_fun(time_marks - a$days[i])
  trigger_constant <- theta[as.numeric(a[i, 'e_type']) + 1]
  output <- colSums(trigger_constant * (g_at_time_marks %o% h_all_locs)) # gives a T x S matrix
  return(output)
}

lambda_at_locs <- function(xs, ys) {
  
  time_integral_bg <- mean(trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks)) * TT
  
  pb <- txtProgressBar(min=1, max=nrow(da), initial=1)
  trig_part_all_times <- Matrix(0, 1, length(xs))
  for (i in 1:nrow(da)) {
    trig_part_all_times <- trig_part_all_times + triggering_at_locs_with_all_times_fun(a = da,
                                                                                       x = xs,
                                                                                       y = ys,
                                                                                       i = i, theta = theta)
    setTxtProgressBar(pb, i)
  }
  
  trigger_contribution <- (TT / length(time_marks)) * trig_part_all_times # This should leave a vector of size S
  lambda_at_probes <- mu0 * background_fun(xs, ys) * time_integral_bg + trigger_contribution
  as.numeric(lambda_at_probes)
}

xs <- c(ripley_bg_basex)[c(as.matrix(ripley_bg_marks)) == TRUE]
ys <- c(ripley_bg_basey)[c(as.matrix(ripley_bg_marks)) == TRUE]

ripley_lambda_all <- lambda_at_locs(xs, ys)
lambda_star <- min(ripley_lambda_all)
harmonic_mean_sq <-  1 / mean(1 / (ripley_lambda_all ^ 2))

total_area <- st_area(shp)

ripley_lambda_obs <- lambda_at_locs(da[da$e_type == 0,]$coorx, da[da$e_type == 0,]$coory)
# Create a distances matrix
orig_events_ids <- which(da$e_type == 0)

n_events <- length(orig_events_ids)
dis_matrix <- matrix(nrow=n_events, ncol=n_events)

for(i in 1:n_events) {
  for(j in 1:n_events) {
    idx_i <- orig_events_ids[i]
    idx_j <- orig_events_ids[j]
    dis_matrix[i, j] <- sqrt((da$coorx[idx_i] - da$coorx[idx_j]) ^ 2 + (da$coory[idx_i] - da$coory[idx_j]) ^ 2)
  }
}

# replicate the weights n_events times to obtain a matrix and zero out the diagonals on it
weights <- lambda_star / as.numeric(ripley_lambda_obs)
weight_mat_repl <- t(replicate(n_events, weights))
diag(weight_mat_repl) <- 0


weighted_ripley <- function(h) {
  # make a copy of the distance matrix and zero-out the far ones
  dis_matrix_close <- dis_matrix
  dis_matrix_close[dis_matrix_close > h] <- 0
  dis_matrix_close[dis_matrix_close > 0] <- 1
  output <- sum(weights * rowSums(dis_matrix_close * weight_mat_repl))
  output <- output / ((lambda_star ^ 2) * total_area)
  # output <- sqrt(output / pi) - h  # variance correction and subtracting the expected value
  output
}

ripley_hs <- seq(0, 0.5, by = 0.0001)
ripley_vals <- numeric(0)
ripley_means <- numeric(0)
ripley_stds <- numeric(0)

for (h in ripley_hs) {
  ripley_vals <- c(ripley_vals, weighted_ripley(h))
  ripley_means <- c(ripley_means, pi * h^2)
  ripley_stds <- c(ripley_stds, sqrt((2 * pi * h^2) / (total_area * harmonic_mean_sq)))
}

ripley_df <- data.frame(x = ripley_hs, y = ripley_vals,
                        mean = ripley_means,
                        std = ripley_stds,
                        high=ripley_means + qnorm(0.975, mean = numeric(length(ripley_hs)), sd=ripley_stds),
                        low=ripley_means + qnorm(0.025, mean = numeric(length(ripley_hs)), sd=ripley_stds))


# Save the data frame if we want to amend the plot
save(ripley_df,file="results/ripley_results.Rda")

p <- ggplot(ripley_df, aes(x, y)) + geom_line() +
  geom_line(aes(x=x, y=mean), linetype = "dotted") + 
  geom_ribbon(aes(ymin = low, ymax = high), alpha=0.1) +
  labs(x = "h [m]", y = parse(text = "K(h)"))
 
ggsave(paste0("figures/", experiment_id, "_ripley_fit_diagnostic.pdf"), plot = p, width = 5.5, height = 3.5)
print(p)