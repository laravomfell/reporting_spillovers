library(ggplot2)
library(scales)
library(viridis)
library(sf)

#text = element_text(size = 14),

theme_print <- function(...){
  theme_classic() +
    theme(strip.text = element_text(colour = 'black'),
          axis.text = element_text(color = "black")) +
    theme(...)
}

theme_set(theme_print())

# plotting mu_bg and overlaying it with streets
mu_bg <- cbind(expand.grid(x = background_base$x,
                          y = background_base$y),
               z = as.vector((background_basevalue * (background_marks)/TT/shp_area)))
# drop undefined values
mu_bg <- mu_bg[mu_bg$z > 0,]


# define a custom labelling function for the scale to turn 0.001 into 10^-3
label_custom <- function(...){
  function(x) label_parse()(gsub("1e", "10^", scientific_format()(x)))
}

p <- ggplot() + 
  # plot underlying background intensity
  geom_tile(data = mu_bg, aes(x * 1000, y * 1000, fill = z)) + 
  # plot colors on log10 scale
  scale_fill_viridis_c(trans = scales::log10_trans(), 
                       labels = label_custom(),
                       direction = -1, option = "magma",
                       breaks = 10^(-c(2,3,4,7)),
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
ggsave("figures/mu_bg.png", plot = p, width = 3.5, height = 4)

# plotting mu_trend
trend <- data.frame(x = time_marks, y = trend_basevalue)

p <- ggplot(trend, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = c(0, 90, 181, 273, 364),
                     labels = c("01/Jan 2018", "01/Apr 2018", "01/Jul 2018", "01/Oct 2018", "31/Dec 2018")) +
  labs(x = "", y = parse(text = "mu[trend](t)")) + theme(panel.grid = element_blank()) +
  theme(plot.margin = margin(t = 4, r = 12, b = 4, l = 4, unit = "pt"))
ggsave("figures/mu_trend.png", plot = p, width = 5.5, height = 3.5)

# plotting mu_weekly
a$datetime_unif <- as.POSIXct(a$datetime_unif, tz = "Europe/London")
wkd_labs <- unique(weekdays(a$datetime_unif, abbreviate = T))
# careful this depends on the data
weekly <- data.frame(x = weekly_base, y = weekly_basevalue)
p <- ggplot(weekly, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = seq(0, 6), labels = wkd_labs, expand = c(0.02, 0.02)) +
  labs(x = "Day of week", y = parse(text = "mu[weekly](t)")) + theme(panel.grid = element_blank())
ggsave("figures/mu_weekly.png", plot = p, width = 5.5, height = 3.5)

# plotting mu_daily
daily <- data.frame(x = daily_base, y = daily_basevalue)
p <- ggplot(daily, aes(x, y)) + geom_line() +
  scale_x_continuous(breaks = seq(0, 1, .25), 
                     labels = c("0:00", "06:00", "12:00", "18:00", "0:00"),
                     expand = c(0.02, 0.02)) +
  theme(panel.grid = element_blank()) +
  labs(x = "Time of day", y = parse(text = "mu[daily](t)"))
ggsave("figures/mu_daily.png", plot = p, width = 5.5, height = 3.5)

# plotting g(t)
gt <- data.frame(x = g_base, y = g_basevalue)
p <- ggplot(gt, aes(x,y)) + geom_point(data = data.frame(mids = g_temp$mids,
                                                    density = g_temp$density), 
                                  aes(mids, density), alpha = .4, shape = 21, color = "#4c4c4c") + 
  geom_line(color = "black", size = 0.9) +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(breaks = c(seq(0, 5, 1), 10, 15)) +
  labs(x = "Time in days", y = "g(t)")
ggsave("figures/gt.png", plot = p, width = 4.5, height=3.5)


# plotting h(s)
h_s <- cbind(expand.grid(x = h_base_x * 1000, y = h_base_y * 1000),
            z = as.vector(h_basevalue))

cont_label <- data.frame(x = -10, 
                         y = c(-180, -140, -110, -80, -35), 
                         label = seq(0.25, 0.65, by = 0.1))

p <- ggplot(h_s[h_s$z > 0.2,], aes(x,y,z = z)) + 
  geom_contour(color = "black", breaks = cont_label$label) +
  labs(x = "Distance in X (in m)", y = "Distance in Y (in m)") + 
  coord_fixed(xlim = c(-200, 200), ylim = c(-200, 200)) + 
  # hacky solution to overplot white on the black circles
  geom_label(data = cont_label, aes(x,y,label = label), 
             inherit.aes = F, 
             color = "white", 
             label.padding = unit(0.15, "lines"),
             size = 2.75) + 
  # plot text on top
  geom_text(data = cont_label, aes(x, y, label = label), inherit.aes = F, size = 2.75) + 
  geom_point(data = data.frame(x = 0, y = 0, z = 0), color = "#4c4c4c", shape = 3)


ggsave("figures/hs.png", plot = p, width = 3.5, height = 3.5)


# multiply through day and week
day_week <- data.frame(x = weekly_base,
                       y = daily_fun(weekly_base) * weekly_fun(weekly_base))

p <- ggplot(day_week, aes(x, y))+ geom_line() +
  scale_x_continuous(breaks = seq(0, 6), labels = paste(wkd_labs, "00:00"), expand = c(0.02, 0.02)) +
  labs(x = "Time", y = parse(text = "mu[daily](t)~x~mu[weekly](t)"))
ggsave("figures/day_week.png", plot = p, width = 5.5, height = 3.5)
