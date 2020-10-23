
 
 bgrates.at.all.locations.no.mu <- ((trend.fun(time.marks)*weekly.fun(time.marks)*daily.fun(time.marks))*
                          mean(background.spatial.fun(background.basex,background.basey)*background.marks)*
                          (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))

triggers.at.all.locations.no.A <- 0

for(i in 1:nrow(a)){ 
           if(as.integer(i/100)*100==i)print(i)
        
           temp <- excite.spatial.fun(background.basex-a$coorx[i], 
                                               background.basey-a$coory[i])
    
           triggers.at.all.locations.no.A <- (triggers.at.all.locations.no.A+ (excite.temporal.fun(time.marks-a$days[i]))*
      mean((background.marks>0)*temp)* (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))
       }

  lambda.at.all.locations <-  mu* bgrates.at.all.locations.no.mu + A * triggers.at.all.locations.no.A


 
postscript("raw-residual.ps", paper="special", height=7, width=14)

 par(mfrow=c(1,2),cex=1.2)

  plot(c(0,a$days), 0:nrow(a), type="s",xlab="Time (days)", ylab="Sum. Freq.", lwd=2, xaxs="i", yaxs="i",xlim=c(0,735), ylim=c(0, 5200),main="(a)",cex.main=1.5)
  abline(0, nrow(a)/730,lwd=1.5, lty=2)   

  plot( cumsum(lambda.at.all.locations)*0.005, stepfun(a$days, (0:nrow(a)))(time.marks), type="s",lwd=2,
       xlab="Transformed time", ylab="Cum. Freq.", xaxs="i", yaxs="i",xlim=c(0,5100), ylim=c(0, 5200),main="(b)",
       cex.main=1.5)
  abline(0,1,lty=2, lwd=1.5)

  points(0:5089, 5089* qbeta(0.975, 0:5089+1, 5089-(0:5089)+1),type="l", lty=2)
  points(0:5089, 5089* qbeta(0.025, 0:5089+1, 5089-(0:5089)+1),type="l", lty=2)

dev.off()
