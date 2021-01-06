library(readxl)
library(data.table)

dens = read_xlsx("d:/research/westmids/data/census/qs102ew_population_density/mid2018.xlsx", sheet = "Mid-2018 Persons", skip = 4)
setDT(dens)
setnames(dens, c("OA11CD", "All Ages"), c("oa_2011", "all_ages"))

lookup = fread("d:/research/westmids/output/census/lookup.csv")
lookup[, area_h := area * 100]

dens <- merge(lookup, dens[, .(oa_2011, all_ages)])
dens[, density_per_hectare := all_ages/area_h]

load("results_trend14/influ_mat.Rdata")
load("results_trend14/theta.Rdata")

a$trig_pres = rowSums(influ_events) * theta[a$e_type + 1]

a_plot = merge(a, dens, by = "oa_2011")


label_10 <- function(...){
  function(x) label_parse()(gsub("e", "~x~10^", gsub("0e\\+00", "0", scientific_format()(x))))
}

p <- ggplot(a_plot, aes(density_per_hectare, trig_pres, color = factor(e_type))) + 
  geom_jitter(alpha = .4) +
  labs(x = "Population density per square hectare", y = "Triggering pressure\nexerted by event") +
  scale_color_manual(values = c("0" = "#721F81FF",
                                "1" = "#E65100"),
                     name = "Type of event",
                     labels = c("Initial report", "Follow-up visit")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(labels = label_10(), breaks = c(0, 2e-07, 4e-07, 6e-07, 8e-07)) +
  guides(color = guide_legend(title.position = "top", override.aes = list(alpha = 1)))
ggsave("figures/dens.png", plot = p, width = 5.5, height = 3.5)


### version 2:

detached <- read_xlsx("ks401ew_dwellings.xlsx", skip = 8)
setDT(detached)
setnames(detached, c("oa_2011", "all", "detached", "semi_detached", "terraced", "flat1", "flat2", "flat3"))
detached <- detached[!is.na(all)]
detached[, share_detached := (detached + semi_detached)/all]

a_plot <- merge(a_plot, detached[, .(oa_2011, share_detached)], by = "oa_2011")

p <- ggplot(a_plot, aes(1 - share_detached, trig_pres)) + 
  geom_point(aes(color = factor(e_type)), alpha = 0.4) + 
  scale_y_continuous(labels = label_10(), breaks = c(0, 2e-07, 4e-07, 6e-07, 8e-07)) +
  scale_color_manual(values = c("0" = "#721F81FF",
                                "1" = "#E65100"),
                     name = "Type of event",
                     labels = c("Initial report", "Follow-up visit")) +
  theme(legend.position = "bottom",
        panel.grid = element_blank()) +
  labs(x = "Share of households living in non-detached houses\nin neighborhood of event",
       y = "Mean number of\nreports triggered") +
  guides(color = guide_legend(title.position = "top",
                              override.aes = list(alpha = 1))) + geom_smooth(method = "lm", color = "black")
ggsave("figures/detached.png", width = 5.5, height = 3.5)
