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
ggsave(paste0("figures/", experiment_id, "_time.pdf"), plot = p, width = 6.5, height = 3.5)