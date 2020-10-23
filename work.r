require('polyCub')
require('spatstat')
require('fields')
source('~/Share/Mylib/Rsrc/inpoly.r')

postscript("functions2.ps", paper="special", height=10, width=10)

# par(mar=c(4,4,1.5,1), mgp=c(1.5,0.5,0), cex.lab=1.25,new=F)

   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
 plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
     title(main="(a)                  ");
 },
 key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2})))} 
)
dev.off()

postscript("functions.ps", paper="special", height=15, width=10)

# x11()

 par(fig=c(0, 1,0.7,1), mar=c(4,4,1.5,1), mgp=c(1.5,0.5,0), cex.lab=1.25,new=F)

 plot(trend.base, trend.basevalue,type="l", lwd=2, col=2, xlab="Time in days",
       ylab=expression(mu[t](t)), main="(b) Trend")


par(fig=c(0, 0.7,0.5,0.7),  new=T)

plot(weekly.base, weekly.basevalue,type="l", lwd=2, col=2, xlab="Time in days", ylab=expression(mu[w](t)), main='(c) Weekly periodicity')



par(fig=c(0.7, 1,0.5,0.7),  new=T)

plot(daily.base, daily.basevalue,type="l", lwd=2, col=2, xlab="Time in days",
      ylab=expression(mu[d](t)),main="(c) Daily peri.")

par(fig=c(0., .475,0.2,0.5),  new=T)

plot(excite.temporal.base, excite.temporal.basevalue,type="l",xlim=c(0,5),xlab="Time in days", ylab=expression(g(t)),col=2, lwd=2, main="(d) Temporal response")
points(excite.temporal.series$mids, excite.temporal.series$density)
 points(excite.temporal.base, excite.temporal.basevalue,type="l",xlim=c(0,5),xlab="Time in days", ylab="Events per day",col=2, lwd=2)

par(fig=c(.525, 1,0.2,0.5),  new=T)

contour(excite.spatial.base.x, excite.spatial.base.y,
     excite.spatial.basevalue, xlab="X (100 m)", ylab="Y (100 m)" ,xlim=c(-1,1), ylim=c(-1,1), main=paste("(e) Spatial response:",expression(h(x,y))))


dev.off()
