# In this file, we plot a time series of the event types and a spatial map of 
# where the events are located
setDT(da)
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
ggsave(paste0("figures/", experiment_id, "_time.pdf"), plot = p, width = 6.5, height = 3.5)

# here, we would normally load our shapefile, instead we just take the domain on 
# which we generated the data

p <- ggplot() +
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
ggsave(paste0("figures/", experiment_id, "_map.pdf"), plot = p, width = 3.5, height = 4) 
