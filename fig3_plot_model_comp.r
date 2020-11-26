library(ggplot2)
library(viridis)
library(osmdata)
library(magrittr)

theme_print <- function(...){
  theme_light() +
    theme(text = element_text(size = 14),
          plot.background = element_blank(),
          legend.background = element_blank(), 
          legend.key = element_blank(),
          strip.text = element_text(colour = 'black'),
          strip.background = element_rect(fill = "lightgray")) +
    theme(...)
}

theme_set(theme_print())

# plotting mu_bg and overlaying it with streets
mu_bg <- cbind(expand.grid(x = background_base$x,
                          y = background_base$y),
               z = as.vector((background_basevalue * (background_marks)/TT/shp_area)))
# drop undefined values
mu_bg <- mu_bg[mu_bg$z > 0,]


if (file.exists("roads.shp")){
  roads <- st_read("roads.shp", crs = 27700)
} else {
  # get streets
  street_values <- c("motorway",
                     "trunk",
                     "primary", "secondary", "tertiary")
  street_values <- c(street_values, paste0(street_values, "_link"),
                     "unclassified", "residential", "living_street", "pedestrian")
  
  
  osm_bbox <- st_bbox(st_transform(st_as_sfc(bbox * 1000), 4326))
  # define query for all roads
  roads <- opq(bbox = osm_bbox, timeout = 50) %>%
    add_osm_feature(key = "highway", value = street_values) %>%
    osmdata_sf()
  roads <- roads$osm_lines
  
  # keep only the lines inside shp
  roads <- st_transform(roads, crs = 27700)
  roads <- st_intersection(roads, shp)
  
  # define line thinkness for streets
  roads$lw <- dplyr::case_when(roads$highway %in% c("residential", 
                                                    "unclassified", 
                                                    "tertiary", 
                                                    "tertiary_link", 
                                                    "living_street") ~ "B",
                               roads$highway == "pedestrian" ~ "C",
                               TRUE ~ "A")
  
  roads <- subset(roads, select = c("osm_id", "name", "highway", "lw"))
  st_write(roads, dsn = "roads.shp", delete_dsn=TRUE)
}


p <- ggplot() + 
  # plot underlying background intensity
  geom_tile(data = mu_bg, aes(x * 1000, y * 1000, fill = z)) + 
  # plot colors on log10 scale
  scale_fill_viridis_c(trans = scales::log10_trans(), direction = -1, option = "magma",
                       breaks = c(0.000001, 0.00001, 0.01),
                       labels = c(0.000001, 0.00001, 0.01)) + 
  # add roads
  geom_sf(data = roads, aes(size = lw), inherit.aes = F) +
  scale_size_manual(values = c("A" = .75, "B" = .65, "C" = .4), guide = FALSE) +
  coord_sf(datum = NA) +
  labs(x = "", y = "")
ggsave("figures/mu_bg.png", plot = p)

# plotting mu_trend
trend <- data.frame(x = time_marks, y = trend_basevalue)

p <- ggplot(trend, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = c(31, 213, 397),
                     labels = c("01/Jan 2016", "01/Jul 2016", "01/Jan 2017")) +
  scale_y_continuous(limits = c(0.4, 1.25)) +
  labs(x = "", y = "mu trend")
ggsave("figures/mu_trend.png", plot = p)

# plotting mu_weekly
# careful this depends on the data
weekly <- data.frame(x = weekly_base, y = weekly_basevalue)
p <- ggplot(weekly, aes(x,y)) + geom_line() +
  scale_x_continuous(breaks = seq(0, 6), labels = c("Tue", "Wed", "Thu", "Fri", "Sat", "Sun", "Mon")) +
  labs(x = "Weekday", y = "mu weekly")
ggsave("figures/mu_weekly.png", plot = p)

# plotting mu_daily
daily <- data.frame(x = daily_base, y = daily_basevalue)
p <- ggplot(daily, aes(x, y)) + geom_line() +
  scale_x_continuous(breaks = seq(0,1, .25), labels = c("0:00", "06:00", "12:00", "18:00", "0:00")) +
  labs(x = "Time of day", y = "mu daily")
ggsave("figures/mu_daily.png", plot = p)

# plotting g(t)
gt <- data.frame(x = g_base, y = g_basevalue)
p <- ggplot(gt, aes(x,y)) + geom_point(data = data.frame(mids = g_temp$mids[-1],
                                                    density = g_temp$density[-1]), 
                                  aes(mids, density), alpha = .4, shape = 21, color = "black") + 
  geom_line(color = "darkgreen", size = 0.9) +
  labs(x = "Time in days", y = "Temporal excitation")
ggsave("figures/gt.png", plot = p)


p <- ggplot(gt, aes(x,y)) + geom_point(data = data.frame(mids = g_temp$mids[-1],
                                                         density = g_temp$density[-1]), 
                                       aes(mids, density), alpha = .4, shape = 21, color = "black") + 
  geom_line(color = "darkgreen", size = 0.9) +
  labs(x = "Time in days", y = "Temporal excitation") +
  coord_cartesian(xlim = c(0, 3))
ggsave("figures/gt_zoomed.png", plot = p)

# plotting h(s)
h_s <- cbind(expand.grid(x = h_base_x * 1000, y = h_base_y * 1000),
            z = as.vector(h_basevalue))

p <- ggplot(h_s, aes(x,y,z =z )) + geom_contour(color = "black") +
  labs(x = "Distance in m", y = "Distance in m")
ggsave("figures/hs.png", plot = p)
