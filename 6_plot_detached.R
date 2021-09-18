# In this file, we create Figure 2.4 where we show the triggering pressure
# against the dominant household type in the area

# first, calculate how many other events this event probably triggered
da$trig_pres <- rowSums(influ_events) * theta[a$e_type + 1]

# data on household type in 
# https://www.nomisweb.co.uk/census/2011/ks401ew
detached <- read_xlsx("ks401ew_dwellings.xlsx", skip = 8)
setDT(detached)
setnames(detached, c("oa_2011", "all", "detached", "semi_detached", "terraced", "flat1", "flat2", "flat3"))
detached <- detached[!is.na(all)]
# calculate share of detached houses
detached[, share_detached := (detached + semi_detached)/all]

# merge in by oa_2011 (not atually possible in our toy data since we didn't create any fake output area labels)
da <- merge(da, detached[, .(oa_2011, share_detached)], by = "oa_2011")

p <- ggplot(da, aes(1 - share_detached, trig_pres)) + 
  geom_point(aes(color = factor(e_type)), alpha = 0.4) + 
  scale_y_continuous(labels = label_10(), breaks = c(0, 5e-09, 1e-08)) +
  scale_color_manual(values = c("0" = "#721F81FF",
                                "1" = "#E65100"),
                     name = "Type of event",
                     labels = c("Initial report", "Follow-up visit")) +
  theme(legend.position = "bottom",
        panel.grid = element_blank()) +
  labs(x = "Share of households living in non-detached houses\nin neighborhood of event",
       y = "Number of\nreports triggered") +
  guides(color = guide_legend(title.position = "top",
                              override.aes = list(alpha = 1))) + 
  geom_smooth(aes(group = factor(e_type)), method = "lm", color = "black")
ggsave(paste0("figures/", experiment_id, "_detached.pdf"), plot = p, width = 5.5, height = 3.5)
