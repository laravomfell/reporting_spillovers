# create figures 2.2 a) through e)

# I want all mu plots to be on the same axis scale, the most variable is day_week so 
# precalculate the ylimits first

# multiply through day and week
day_week <- data.frame(x = weekly_base,
                       y = daily_fun(weekly_base) * weekly_fun(weekly_base))
# get ylims
y_lim <- range(day_week$y)
y_breaks <- pretty(y_lim)

# a) mu_trend
trend <- data.frame(x = time_marks, y = trend_basevalue)

p <- ggplot(trend, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = c(0, 90, 181, 273, 364),
                     labels = c("01/Jan 2018", "01/Apr 2018", "01/Jul 2018", 
                                "01/Oct 2018", "31/Dec 2018")) +
  labs(x = "", y = parse(text = "mu[trend](t)")) + 
  scale_y_continuous(limits = y_lim, breaks = y_breaks) +
  theme(plot.margin = margin(t = 4, r = 12, b = 4, l = 4, unit = "pt"))
ggsave(paste0("figures/", experiment_id, "_mu_trend.pdf"), plot = p, width = 5.5, height = 3.5)


# b) mu_daily
daily <- data.frame(x = daily_base, y = daily_basevalue)
p <- ggplot(daily, aes(x, y)) + geom_line() +
  scale_x_continuous(breaks = seq(0, 1, .25), 
                     labels = c("00:00", "06:00", "12:00", "18:00", "00:00"),
                     expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = y_lim, breaks = y_breaks) +
  labs(x = "Time of day", y = parse(text = "mu[daily](t)"))
ggsave(paste0("figures/", experiment_id, "_mu_daily.pdf"), plot = p, width = 5.5, height = 3.5)


# c) mu_weekly

da$datetime_unif <- as.POSIXct(da$datetime_unif, tz = "Europe/London")
wkd_labs <- unique(weekdays(da$datetime_unif, abbreviate = T))

weekly <- data.frame(x = weekly_base, y = weekly_basevalue)
p <- ggplot(weekly, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = seq(0, 6), labels = wkd_labs, expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = y_lim, breaks = y_breaks) +
  labs(x = "Day of week", y = parse(text = "mu[weekly](t)"))
ggsave(paste0("figures/", experiment_id, "_mu_weekly.pdf"), plot = p, width = 5.5, height = 3.5)



# d) day-week
p <- ggplot(day_week, aes(x, y))+ geom_line() +
  scale_x_continuous(breaks = seq(0, 6), 
                     labels = paste(wkd_labs, "00:00"), 
                     expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = y_lim, breaks = y_breaks) +
  labs(x = "Time", y = parse(text = "mu[daily](t)~x~mu[weekly](t)"))
ggsave(paste0("figures/", experiment_id, "_day_week.pdf"), plot = p, width = 5.5, height = 3.5)


# e) mu_area
mu_area <- cbind(expand.grid(x = background_base$x,
                             y = background_base$y),
                 z = as.vector((background_basevalue * (background_marks))))
# drop undefined values
mu_area <- mu_area[mu_area$z > 0,]

p <- ggplot() + 
  # plot underlying background intensity
  geom_tile(data = mu_area, aes(x, y, fill = z)) + 
  # plot colors on log10 scale
  scale_fill_viridis_c(trans = scales::log10_trans(), 
                       labels = label_custom(),
                       direction = -1, option = "magma",
                       #breaks = c(0.01, 1, 100), # TODO: removed temporarily
                       name = parse(text = "mu[area](s)")) + 
  coord_sf(datum = NA, expand = FALSE) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", 
        panel.border = element_blank(),
        legend.margin=margin(t = -.75, unit='cm'),
        plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit="pt"),
        legend.title=element_text(size=12),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(title.position = "top",
                                barheight = unit(4, units = "mm"),
                                barwidth = unit(45, units = "mm")))
ggsave(paste0("figures/", experiment_id, "_mu_bg.png"), plot = p, width = 3.5, height = 4)

# 2.3

# a) g(t)
gt <- data.frame(x = g_base, y = g_basevalue)
p <- ggplot(gt, aes(x,y)) + geom_point(data = data.frame(mids = g_temp$mids,
                                                    density = g_temp$density), 
                                  aes(mids, density), alpha = .4, shape = 21, color = "#4c4c4c") + 
  geom_line(color = "black", size = 0.9) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(seq(0, 5, 1), 10, 15, 20, 25, 30)) +
  labs(x = "Time in days", y = "g(t)") + geom_hline(yintercept = 0, linetype = "dotted")
ggsave(paste0("figures/", experiment_id, "_gt.pdf"), plot = p, width = 4.5, height=3.5)


# b) h(s)
h_s <- cbind(expand.grid(x = h_base_x * 1000, y = h_base_y * 1000),
             z = as.vector(h_basevalue))
## cont_label <- data.frame(x = -10, 
##                          y = c(-315, -252, -192,-125, 0), 
##                          b = seq(0.12, 0.16, by = 0.01),
##                          label = c(seq(0.12, 0.15, by = 0.01), ""))

## p <- ggplot(h_s[h_s$z >= 0.1,], aes(x,y,z = z)) + 
##   geom_contour(color = "black", breaks = cont_label$b) +
##   labs(x = "Distance in X (in m)", y = "Distance in Y (in m)") + 
##   coord_fixed(xlim = c(-500, 350), ylim = c(-400, 400)) + 
##   # hacky solution to overplot white on the black circles
##   geom_label(data = cont_label, aes(x,y,label = label), 
##              inherit.aes = F, 
##              color = "white", 
##              label.padding = unit(0.15, "lines"),
##              size = 2.75) + 
##   # plot text on top
##   geom_text(data = cont_label, aes(x, y, label = label), inherit.aes = F, size = 2.75) + 
##   geom_point(data = data.frame(x = 0, y = 0, z = 0), color = "#4c4c4c", shape = 3)

# This is temporary for the experimentation stage.
p <- ggplot(h_s, aes(x, y, z= z, colour=stat(level))) +
    labs(x = "Distance in X (in m)", y = "Distance in Y (in m)") + 
    geom_contour() +
    scale_color_gradient(low = "gray100", high = "gray0", space = "Lab", 
                         guide = guide_colorbar(title="h(s)", title.position="top"))                         
ggsave(paste0("figures/", experiment_id, "_hs.pdf"), plot = p, width = 3.5, height = 3.5)
