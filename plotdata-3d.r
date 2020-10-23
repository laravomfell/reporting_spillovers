rgl.close()
require('rgl')

plot3d(a$coorx, a$coory, -a$days, xlab="", ylab="", zlab="",col=rainbow(760)[a$days],xlim=c(749.3, 754.4), ylim=c(4428.5, 4432.6), zlim=c(-730,0), type="s",
       radius=3, axes=F)

axes3d(font=3)
box3d()
#axis3d('x', pos = c(NA, min(pretty(a$coory)), 0))
#axis3d('x', pos = c(NA, max(pretty(a$coory)), 0))


#axis3d("x", pos=c(NA,0,0))
title3d('', '', 'X (100m)', 'Y (100m)', 'Z (days)')
   
#       rgl.bringtotop()
  
lines3d(city.boundary$x, city.boundary$y, city.boundary$x*0-0, lty=1,lwd=3)

lines3d(poly.downtown$V1, poly.downtown$V2, poly.downtown$V1*0-0, lty=1,col=1,lwd=3)
