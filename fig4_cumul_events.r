no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

bg_at_all_locations <- trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks) * mean(background_fun(background_basex, 
                                                                                                                    background_basey)[as.vector(background_marks > 0)]) * ra * bg_weight


trigger_at_all_locations <- foreach(i = 1:nrow(a)) %dopar% trigger_at_all_fun(i = i, constants = bg_weight * ra)
# reduce by event type
trigger_at_all_locations <- map(event_types, function(x) reduce(trigger_at_all_locations[a$e_type == x],
                                                                `+`))

# multiply through with theta and simplify
lambda_at_all_locations <- mu0 * bg_at_all_locations + reduce(map2(trigger_at_all_locations, theta, `*`), `+`)
save(lambda_at_all_locations, file = "results/lambda_at_all_locations.Rdata")

# creating fig4(b)
n_events <- nrow(a[a$e_type == 0,])

fit <- data.frame(y = cumsum(lambda_at_all_locations) * space_units,
                  x = stepfun(a$days[a$e_type == 0], 0:n_events)(time_marks))
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
  scale_x_continuous(breaks = c(0, 2000, 4000, n_events)) + 
  scale_y_continuous(breaks = c(0, 2000, 4000, n_events)) +
  coord_cartesian(xlim = c(0, max(fit$x)), ylim = c(0, n_events))+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9)) +
  labs(x = "Event index", y = "Cumulative frequency")
ggsave("figures/trans_time.png", plot = p, width = 5.5, height = 3.5)

# trans time 2
conf <- data.frame(x = 0:n_events, y = 0,
                   ymax = n_events * qbeta(.975, 0:n_events + 1, n_events - (0:n_events) + 1),
                   ymin = n_events * qbeta(.025, 0:n_events + 1, n_events - (0:n_events) + 1))

fit$dev <- fit$y - fit$x

p <- ggplot(fit, aes(y, dev)) + 
  geom_ribbon(data = conf,
              aes(x, ymin = ymin -x, ymax = ymax -x, y=y),inherit.aes = F, fill = "#d6d6d6") +
  geom_hline(yintercept = 0) + 
  scale_x_continuous(breaks = c(0, 2000, 4000, n_events)) + 
  scale_y_continuous(breaks = c(-50,-25,0,25,50)) +
  geom_line(color = "#721F81FF") +
  labs(x = "Event index", y = parse(text = "Predicted~-~expected~number~of~events"))
ggsave("figures/dev.png", plot = p, width = 5.5, height = 3.5)
##### trans resid

tri_i = cumsum(lambda_at_all_locations) * space_units
tau_i = tri_i - shift(tri_i)
z_i = 1 - exp(-tau_i)
zi_sort = sort(z_i)
plot((1:36500 - 0.5)/36500, zi_sort)

# hm.
