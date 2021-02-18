# In this file, we create Figure 2.2f) which shows domestic abuse against deprivation
# again, this code cannot actually be used since we did not create any fake lower layer super output area labels


# calculate background area intensity at event
da[e_type == 0, bg_area := background_fun(coorx, coory)]

# deprivation data available from 
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/467764/File_1_ID_2015_Index_of_Multiple_Deprivation.xlsx
deprivation <- read_xlsx("File_1_ID_2015_Index_of_Multiple_Deprivation.xlsx", sheet = "IMD 2015")
setnames(deprivation, snakecase::to_any_case(colnames(deprivation)))
da <- merge(da, deprivation, by.x = "lsoa_2011", by.y = "lsoa_code_2011")

# need to manually add jitter to spread out points more
da[, mu_bg := 10^(log10(bg_area) + runif(.N, -0.02, 0.02))]
da[, depr := multiple_depr_score + runif(.N, -0.6, 0.6)]

p <- ggplot(da[e_type == 0], 
            aes(mu_bg, index_of_multiple_deprivation_imd_rank_where_1_is_most_deprived, 
                color = mu_bg)) + 
  geom_point(alpha = .4) +
  geom_point(shape = 21) +
  scale_x_continuous(trans = scales::log10_trans(), 
                     labels = label_custom(),
                     breaks = c(0.01, 1, 100)) +
  # plot colors on log10 scale
  scale_color_viridis_c(trans = scales::log10_trans(), 
                        direction = -1, option = "magma",
                        guide = F) + 
  labs(x = parse(text = "mu[area](s)"),
       y = "Multiple deprivation score")
ggsave("figures/deprivation.pdf", plot = p, width = 5.5, height = 4)
