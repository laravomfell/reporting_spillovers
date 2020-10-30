
source('~/Share/Mylib/Rsrc/inpoly.r')

poly.downtown <- read.table('castellon-city-center.polygon.dat')/1000

postscript('r2.ps', height=8, width=10)

############ Reconstructed temporal response function for 2012 and 2013 #################

if(1==1){
ids <-  a$days[ij.mat[,2]]<= 365

excite.temporal.series.hh <- hist.weighted((a$days[ij.mat[,1]]-a$days[ij.mat[,2]]),
                                           excite.wghs
             /(excite.spatial.edge.correction[ij.mat[,2]]
             *temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue.hh<- ker.smooth.conv(excite.temporal.series.hh$mids, excite.temporal.series.hh$density, bandwidth=0.03)

excite.temproal.basevalue.hh <- excite.temporal.basevalue.hh / simpson(excite.temporal.basevalue, 0.005)

   excite.temporal.series.h1 <- hist.weighted((a$days[ij.mat[,1]]-a$days[ij.mat[,2]]) [ids],
                                           excite.wghs[ids]
             /(excite.spatial.edge.correction[ij.mat[,2]]
             *temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))[ids]
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue.h1<- ker.smooth.conv(excite.temporal.series.h1$mids, excite.temporal.series.h1$density, bandwidth=0.03)

excite.temproal.basevalue.h1 <- excite.temporal.basevalue.h1 / simpson(excite.temporal.basevalue, 0.005)

plot(excite.temporal.series$mids, excite.temporal.series$density, pch=".",cex=2,col=1, main=paste("mu=",mu,"A=",A), xlim=c(0, 2))

points(excite.temporal.base, excite.temporal.basevalue.hh,type="l", lwd=2)


points(excite.temporal.series.h1$mids, excite.temporal.series.h1$density, pch=".",cex=2,col=2, main=paste("mu=",mu,"A=",A), xlim=c(0, 2))

points(excite.temporal.base, excite.temporal.basevalue.h1,type="l",col=2, lwd=2)



  
 excite.temporal.series.h2 <- hist.weighted((a$days[ij.mat[,1]]-a$days[ij.mat[,2]]) [!ids],
                                           excite.wghs[!ids]
             /(excite.spatial.edge.correction[ij.mat[,2]]
             *temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))[!ids]
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue.h2<- ker.smooth.conv(excite.temporal.series.h2$mids, excite.temporal.series.h2$density, bandwidth=0.03)

excite.temproal.basevalue.h2 <- excite.temporal.basevalue.h2 / simpson(excite.temporal.basevalue, 0.005)

points(excite.temporal.base, excite.temporal.basevalue.h2,type="l", col=3,lwd=2)

points(excite.temporal.series.h2$mids, excite.temporal.series.h2$density, pch=".",cex=2, main=paste("mu=",mu,"A=",A), xlim=c(0, 2),col=3)



############ Reconstructed temporal response function for downtown and suburb #################


ids = inpoly( a$coorx, a$coory, poly.downtown$V1, poly.downtown$V2)>=0

   excite.temporal.series.h1 <- hist.weighted((a$days[ij.mat[,1]]-a$days[ij.mat[,2]]) [ids],
                                           excite.wghs[ids]
             /(excite.spatial.edge.correction[ij.mat[,2]]
             *temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))[ids]
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue.h1<- ker.smooth.conv(excite.temporal.series.h1$mids, excite.temporal.series.h1$density, bandwidth=0.03)

excite.temproal.basevalue.h1 <- excite.temporal.basevalue.h1 / simpson(excite.temporal.basevalue, 0.005)

plot(excite.temporal.series$mids, excite.temporal.series$density, pch=".",cex=2,col=1, main=paste("mu=",mu,"A=",A), xlim=c(0, 2))

points(excite.temporal.base, excite.temporal.basevalue.hh,type="l", lwd=2)


points(excite.temporal.series.h1$mids, excite.temporal.series.h1$density, pch=".",cex=2,col=2, main=paste("mu=",mu,"A=",A), xlim=c(0, 2))

points(excite.temporal.base, excite.temporal.basevalue.h1,type="l",col=2, lwd=2)

  
 excite.temporal.series.h2 <- hist.weighted((a$days[ij.mat[,1]]-a$days[ij.mat[,2]]) [!ids],
                                           excite.wghs[!ids]
             /(excite.spatial.edge.correction[ij.mat[,2]]
             *temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))[!ids]
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue.h2<- ker.smooth.conv(excite.temporal.series.h2$mids, excite.temporal.series.h2$density, bandwidth=0.03)

excite.temproal.basevalue.h2 <- excite.temporal.basevalue.h2 / simpson(excite.temporal.basevalue, 0.005)

points(excite.temporal.base, excite.temporal.basevalue.h2,type="l", col=3,lwd=2)

points(excite.temporal.series.h2$mids, excite.temporal.series.h2$density, pch=".",cex=2, main=paste("mu=",mu,"A=",A), xlim=c(0, 2),col=3)

}

###############################Constructing spatial response  for 2012 and 2013 ###############################
ids <-  a$days[ij.mat[,2]]<= 365
excite.spatial.series.h1 <-  hist.weighted.2D(dis.mat[ids,1], dis.mat[ids,2], excite.wghs[ids]
          #/(excite.temporal.edge.correction[ij.mat[,2]])
                                  , 
            x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)
             
excite.spatial.basevalue.h1 <- ker.smooth.2D.fft(excite.spatial.series.h1$x.mids,excite.spatial.series.h1$y.mids, 
 excite.spatial.series.h1$density, x.bandwidth=0.05, y.bandwidth=0.05)

excite.spatial.basevalue.h1 <- excite.spatial.basevalue.h1/simpson.2D(excite.spatial.basevalue.h1,0.002,0.002)

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue.h1[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed',
               xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
    
excite.spatial.series.h2 <-  hist.weighted.2D(dis.mat[!ids,1], dis.mat[!ids,2], excite.wghs[!ids]
          #/(excite.temporal.edge.correction[ij.mat[,2]])
                                  , 
            x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)
             
excite.spatial.basevalue.h2 <- ker.smooth.2D.fft(excite.spatial.series.h2$x.mids,excite.spatial.series.h2$y.mids, 
 excite.spatial.series.h2$density, x.bandwidth=0.05, y.bandwidth=0.05)

excite.spatial.basevalue.h2 <- excite.spatial.basevalue.h2/simpson.2D(excite.spatial.basevalue.h2,0.002,0.002)
 
filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue.h2[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
    

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               (excite.spatial.basevalue.h2-excite.spatial.basevalue.h1)[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), zlim=c(-1,1.5))
    

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               ((excite.spatial.basevalue.h2-excite.spatial.basevalue.h1)/(excite.spatial.basevalue.h2+excite.spatial.basevalue.h1))  [seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), zlim=c(-1,1))
    


###############################Constructing spatial response  for downtown and suburb ###############################

ids = inpoly( a$coorx, a$coory, poly.downtown$V1, poly.downtown$V2)>=0

excite.spatial.series.h1 <-  hist.weighted.2D(dis.mat[ids,1], dis.mat[ids,2], excite.wghs[ids]
          #/(excite.temporal.edge.correction[ij.mat[,2]])
                                  , 
            x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)
             
excite.spatial.basevalue.h1 <- ker.smooth.2D.fft(excite.spatial.series.h1$x.mids,excite.spatial.series.h1$y.mids, 
 excite.spatial.series.h1$density, x.bandwidth=0.05, y.bandwidth=0.05)

excite.spatial.basevalue.h1 <- excite.spatial.basevalue.h1/simpson.2D(excite.spatial.basevalue.h1,0.002,0.002)

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue.h1[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed',
               xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
    
excite.spatial.series.h2 <-  hist.weighted.2D(dis.mat[!ids,1], dis.mat[!ids,2], excite.wghs[!ids]
          #/(excite.temporal.edge.correction[ij.mat[,2]])
                                  , 
            x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)
             
excite.spatial.basevalue.h2 <- ker.smooth.2D.fft(excite.spatial.series.h2$x.mids,excite.spatial.series.h2$y.mids, 
 excite.spatial.series.h2$density, x.bandwidth=0.05, y.bandwidth=0.05)

excite.spatial.basevalue.h2 <- excite.spatial.basevalue.h2/simpson.2D(excite.spatial.basevalue.h2,0.002,0.002)
 
filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue.h2[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5))
    

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               (excite.spatial.basevalue.h2-excite.spatial.basevalue.h1)[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), zlim=c(-2.5,2.5))
    

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               ((excite.spatial.basevalue.h2-excite.spatial.basevalue.h1)/(excite.spatial.basevalue.h2+excite.spatial.basevalue.h1))  [seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed', xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), zlim=c(-1,1))
    


dev.off()
