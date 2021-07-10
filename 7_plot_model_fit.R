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
p

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
p

ggsave(paste0("figures/", experiment_id, "_dev.pdf"), plot = p, width = 5.5, height = 3.5)
