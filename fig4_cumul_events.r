no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

bg_at_all_locations <- trend_fun(time_marks) * weekly_fun(time_marks) * daily_fun(time_marks) * mean(background_fun(background_basex, 
                                                                                                                    background_basey) * background_marks) * ra


trigger_at_all_locations <- foreach(i = 1:nrow(a)) %dopar% trigger_at_all_fun(i = i, constants = bg_weight * ra)
# reduce by event type
trigger_at_all_locations <- map(event_types, function(x) reduce(trigger_at_all_locations[a$e_type == x],
                                                                `+`))

# multiply through with theta and simplify
lambda_at_all_locations <- mu0 * bg_at_all_locations + reduce(map2(trigger_at_all_locations, theta, `*`), `+`)
save(lambda_at_all_locations, file = "results/lambda_at_all_locations.Rdata")


# creating fig 4(a)
p <- ggplot(data.frame(x = c(0, a$days),
                  y = 0:nrow(a)), aes(x,y)) +
  geom_abline(intercept = 0, slope = nrow(a)/TT) +
  geom_line(size = .8, color = "darkgreen")
ggsave("figures/org_time.png", plot = p)

# creating fig4(b)
p <- ggplot(data.frame(x = cumsum(lambda_at_all_locations) * space_units,
                  y = stepfun(a$days, 0:nrow(a))(time_marks)),
       aes(x,y)) + 
  geom_line(color = "darkgreen") + geom_abline(intercept = 0, slope = 1,color = "black") +
  geom_ribbon(data = data.frame(x = 0:nrow(a),y = 0,
                         ymax = nrow(a) * qbeta(.975, 0:nrow(a) + 1, nrow(a) - (0:nrow(a)) + 1),
                         ymin = nrow(a) * qbeta(.025, 0:nrow(a) + 1, nrow(a) - (0:nrow(a)) + 1)),
              aes(x, ymin = ymin, ymax = ymax), alpha = .4) +
  scale_x_continuous(breaks = seq(0, 10000, 2000)) + 
  scale_y_continuous(breaks = seq(0, 10000, 2000))
ggsave("figures/trans_time.png", plot = p)
