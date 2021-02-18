# In this file, we plot a time series of the event types and a spatial map of where the events are located

# drop time stamp
da[, date_unif := as.Date(datetime_unif)]
# turn e_type into a factor so the count operation below includes 0 counts
da[, e_type := factor(e_type, levels = c(0, 1))]
e_count <- da[, .(count = as.vector(table(e_type)), 
                  e_type = levels(da$e_type)), 
              by = date_unif]

p <- ggplot(e_count,
       aes(date_unif, count, color = e_type)) + geom_line(size = .7) +
  labs(x = "", y = "Daily number of events") +
  scale_color_manual(values = c("0" = "#721F81FF",
                                  "1" = "#E65100"),
                     name = "Type of event",
                     labels = c("Initial report", "Follow-up visit")) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.margin=margin(t = -0.5, unit='cm'),
        plot.margin=margin(t = 5.5, r = 5.5, b = 0, l = 5.5, unit="pt")) +
  guides(colour = guide_legend(title.position = "top",
                               override.aes = list(size = 1)))
ggsave("figures/time.pdf", plot = p, width = 6.5, height = 3.5)

# here, we would normally load our shapefile, instead we just take the square on which we generated the data

# our plot also shows the underlying road network if you have a bbox and osmdata
# I uncommented these steps here

# # get streets
# street_values <- c("motorway",
#                    "trunk",
#                    "primary", "secondary", "tertiary")
# street_values <- c(street_values, paste0(street_values, "_link"),
#                    "unclassified", "residential", "living_street", "pedestrian")
#   
#   
# osm_bbox <- st_bbox(st_transform(st_as_sfc(bbox * 1000), 4326))
# # define query for all roads
# roads <- opq(bbox = osm_bbox, timeout = 50) %>%
#   add_osm_feature(key = "highway", value = street_values) %>%
#   osmdata_sf()
# roads <- roads$osm_lines
# 
# # keep only the lines inside shp
# roads <- st_transform(roads, crs = 27700)
# roads <- st_intersection(roads, shp)
# 
# # define line thinkness for streets
# roads$lw <- dplyr::case_when(roads$highway %in% c("residential", 
#                                                   "unclassified", 
#                                                   "tertiary", 
#                                                   "tertiary_link", 
#                                                   "living_street") ~ "B",
#                              roads$highway == "pedestrian" ~ "C",
#                              TRUE ~ "A")
# 
# roads <- subset(roads, select = c("osm_id", "name", "highway", "lw"))


p <- ggplot() +
  # geom_sf(data = roads, aes(size = lw), inherit.aes = F, color = "#4c4c4c") +
  # scale_size_manual(values = c("A" = .65, "B" = .55, "C" = .3), guide = FALSE) +
  geom_sf(data = shp, fill = "transparent") +
  geom_point(data = da, aes(coorx, coory, color = days), alpha = 0.4) +
  scale_color_gradient2(high = "#8a3000", low = "#f5b999", mid = "#e65100", midpoint = 181,
                       name = "Time of event",
                       breaks = c(0, 182, 365), limits = c(0, 365),
                       labels = c("Jan", "Jul", "Dec")) +
  labs(x = "", y = "") +
  coord_sf(datum = NULL, expand = FALSE) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.margin=margin(t = -0.75, unit='cm'),
        plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit="pt"),
        panel.border = element_blank(),
        legend.title=element_text(size=11),
        legend.text = element_text(size = 10),
        axis.line = element_blank()) +
  guides(color = guide_colourbar(title.position = "top",
                                barheight = unit(4, units = "mm"),
                                barwidth = unit(45, units = "mm")))
ggsave("figures/map.pdf", plot = p, width = 3.5, height = 4) 
