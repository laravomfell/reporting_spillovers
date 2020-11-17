library(ggplot2)


# plotting mu_bg

df_mu <- list(x = background_base$x[seq(1, length(background_base$x), 18)],
              y = background_base$y[seq(1, length(background_base$y), 18)],
              z = as.matrix(log10(background_basevalue * (background_marks)/TT/boundary_area))[seq(1, length(background_base$x), 18),
                                                                                      seq(1, length(background_base$y), 18)])

dimnames(df_mu$z) <- list(x = df_mu$x, y = df_mu$y)

m2 = reshape2::melt(df_mu, value.name = "z")
m2 = m2[!is.na(m2$x),]
ggplot(m2, aes(x,y,z = z, fill = z)) + geom_tile() # or geom_raster(interpolate = TRUE)

# plotting mu_trend
# I think this needs more
trend <- data.frame(x = trend_base, y = trend_basevalue)

ggplot(trend, aes(x,y)) + geom_line()
 par(fig=c(0, 1,0.7,1), mar=c(4,4,1.5,1), mgp=c(1.5,0.5,0), cex.lab=1.25,new=F)

 plot(trend.base, trend.basevalue,type="l", lwd=2, col=2, xlab="Time in days",
       ylab=expression(mu[t](t)), main="(b) Trend")


# plotting mu_weekly
weekly <- data.frame(x = weekly_base, y = weekly_basevalue)
ggplot(weekly, aes(x,y)) + geom_line()

# plotting mu_daily
daily <- data.frame(x = daily_base, y = daily_basevalue)
ggplot(daily, aes(x, y)) + geom_line()

# plotting g(t)
gt <- data.frame(x = g_base, y = g_basevalue)
ggplot(gt, aes(x,y)) + geom_points(data = g_temp, aes(mids, density)) + geom_line()

# plotting h(s)
df_hs <- list(x = h_base_x,
              y = h_base_y,
              z = h_basevalue)
dimnames(df_hs$z) = list(x = h_base_x, y = h_base_y)

df_hs <- reshape2::melt(df_hs, value.name = "z")
df_hs <- df_hs[!is.na(df_hs$x),]

ggplot(df_hs, aes(x,y,z =z )) + geom_contour()

