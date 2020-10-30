#(f) Construct spatial component of excitations. Remember 
#\texttt{excite.spatial.base.x}, \texttt{excite.spatial.base.y}, and \texttt{excite.spatial.basevalue}.

dis.mat <- cbind(a$coorx[ij.mat[,1]] - a$coorx[ij.mat[,2]], 
                 a$coory[ij.mat[,1]] - a$coory[ij.mat[,2]])


excite.spatial.series <-  hist.weighted.2D(dis.mat[,1], dis.mat[,2], excite.wghs/(excite.temporal.edge.correction[ij.mat[,2]]), x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)

excite.spatial.mark2 <- (excite.spatial.basex ^2 + excite.spatial.basey^2 < 2.35^2)

temp<- spatial.repetance.fun(excite.spatial.series$x.mids%o%rep(1, length(excite.spatial.series$y.mids)),
                            rep(1, length(excite.spatial.series$y.mids))%o%excite.spatial.series$x.mids)
                             
                         
#excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,excite.spatial.series$y.mids, 
# excite.spatial.series$density/temp, x.bandwidth=0.1, y.bandwidth=0.1)*excite.spatial.mark2

excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,excite.spatial.series$y.mids, 
 excite.spatial.series$density, x.bandwidth=0.1, y.bandwidth=0.1)*excite.spatial.mark2

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue,0.002,0.002)

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed')
