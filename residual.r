
source('~/Share/Mylib/Rsrc/inpoly.r')

poly.downtown <- read.table('castellon-city-center.polygon.dat')/1000

 postscript('residual1.ps',paper='special', height=8, width=8)
# postscript('temp.ps',paper='special', height=8, width=8)
print(par()$mgp)
par(cex=1.5, cex.axis=0.8, mgp=c(1.5,0.5,0))

if (1==1){
#### Detect daily period difference within in the two years ########
new.marks.h1 <-(a$days - as.integer(a$days)) [a$days<=365]
temp <- hist.weighted (new.marks.h1, wghs.daily[a$days<=365], breaks= daily.base)
daily.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h1 <- daily.basevalue.h1/ mean(daily.basevalue.h1)

new.marks.h2 <-(a$days - as.integer(a$days))  [a$days>=365]
temp <- hist.weighted (new.marks.h2, wghs.daily[a$days>=365], breaks= daily.base)
daily.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2 <- daily.basevalue.h2/ mean(daily.basevalue.h2)

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",
     ylab="Relative Level")

points(daily.base, daily.basevalue.h1,type="l",col=2)
points(daily.base, daily.basevalue.h2,type="l",col=3)

 legend(0.15, 1.8, lwd=2, col=1:3,
        legend=c("All", "2012", "2013"),cex=0.8
        )

##### Detect daily periodicity difference in the weekdays ##########

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2),lwd=2, xlab="Time (days)",
     ylab="Relative Level")

for(iwd in 0:6){
   ids <- (as.integer(a$days) - 7*as.integer(a$days/7) == iwd)
 new.marks.hh <- (a$days - as.integer(a$days)) [ids]
temp <- hist.weighted (new.marks.hh, wghs.daily[ids], breaks= daily.base)
daily.basevalue.hh <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.hh <- daily.basevalue.hh/ mean(daily.basevalue.hh)

 points(daily.base, daily.basevalue.hh,type="l",col=iwd+2,lwd=2)
}
 legend(0.13, 1.9, lwd=2, col=c(1,2+0:6),
        legend=c("All", "Su", "Mo", "Tu",                                                    "We", "Th", "Fr", "Sa"), cex=0.7
        )


##### Detect daily periodicity difference in four seasons  ##########

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2),lwd=2, xlab="Time (days)",
     ylab="Relative Level")

for(iwd in 0:3){
   ids <- ( (365/4*iwd < a$days & 365/4* (iwd+1) < a$days) | (365/4*(iwd+4) < a$days & 365/4* (iwd+5) < a$days))
 new.marks.hh <- (a$days - as.integer(a$days)) [ids]
temp <- hist.weighted (new.marks.hh, wghs.daily[ids], breaks= daily.base)
daily.basevalue.hh <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.hh <- daily.basevalue.hh/ mean(daily.basevalue)

 points(daily.base, daily.basevalue.hh,type="l",col=iwd+2,lwd=2)
}

 legend(0.15, 1.8, lwd=2, col=c(1, 2:5),
        legend=c("All", "Jan-Mar", "Apr-Jun", "Jul-Sep",                                                    "Oct-Dec"),cex=0.8
        )


##### Detect daily periodicity difference in space : West & East  ##########

ids = a$coorx <= median(a$coorx)

new.marks.h1 <-(a$days - as.integer(a$days)) [ids]
temp <- hist.weighted (new.marks.h1, wghs.daily[ids], breaks= daily.base)
daily.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h1 <- daily.basevalue.h1/ mean(daily.basevalue.h1)

new.marks.h2 <-(a$days - as.integer(a$days))  [!ids]
temp <- hist.weighted (new.marks.h2, wghs.daily[!ids], breaks= daily.base)
daily.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2 <- daily.basevalue.h2/ mean(daily.basevalue.h2)

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",
     ylab="Relative Level")

points(daily.base, daily.basevalue.h1,type="l",col=2)
points(daily.base, daily.basevalue.h2,type="l",col=3)

 legend(0.15, 1.8, lwd=2, col=1:3,
        legend=c("All", "West", "East"),cex=0.8
        )

##### Detect daily periodicity difference in space : North & Source  ##########

ids = a$coory <= median(a$coory)

new.marks.h1 <-(a$days - as.integer(a$days)) [ids]
temp <- hist.weighted (new.marks.h1, wghs.daily[ids], breaks= daily.base)
daily.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h1 <- daily.basevalue.h1/ mean(daily.basevalue.h1)

new.marks.h2 <-(a$days - as.integer(a$days))  [!ids]
temp <- hist.weighted (new.marks.h2, wghs.daily[!ids], breaks= daily.base)
daily.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2 <- daily.basevalue.h2/ mean(daily.basevalue.h2)

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",
     ylab="Relative Level")

points(daily.base, daily.basevalue.h1,type="l",col=2)
points(daily.base, daily.basevalue.h2,type="l",col=3)

 legend(0.15, 1.8, lwd=2, col=1:3,
        legend=c("All", "South", "North"), cex=.8
        )

##### Detect daily periodicity difference in space : downtown & Suburb  ##########

ids = inpoly( a$coorx, a$coory, poly.downtown$V1, poly.downtown$V2)>=0

new.marks.h1 <-(a$days - as.integer(a$days)) [ids]
temp <- hist.weighted (new.marks.h1, wghs.daily[ids], breaks= daily.base)
daily.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h1 <- daily.basevalue.h1/ mean(daily.basevalue.h1)

new.marks.h2 <-(a$days - as.integer(a$days))  [!ids]
temp <- hist.weighted (new.marks.h2, wghs.daily[!ids], breaks= daily.base)
daily.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2 <- daily.basevalue.h2/ mean(daily.basevalue.h2)

plot(daily.base, daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",
     ylab="Relative Level")

points(daily.base, daily.basevalue.h1,type="l",col=2)
points(daily.base, daily.basevalue.h2,type="l",col=3)

legend(0.1, 1.8, lwd=2, col=c(1:3),
        legend=c("All", "Downtown", "Suburb"),cex=0.7
        )
    
plot(daily.base, daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",
     ylab="Relative Level")


points(daily.base, daily.basevalue.h1,type="l",col=2)
points(daily.base, daily.basevalue.h2,type="l",col=3)

ids2 <- a$coorx < median(a$coorx)
    
new.marks.h2w <-(a$days - as.integer(a$days))  [!ids & ids2]
temp <- hist.weighted (new.marks.h2w, wghs.daily[!ids& ids2], breaks= daily.base)
daily.basevalue.h2w <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2w <- daily.basevalue.h2w/ mean(daily.basevalue.h2w)
points(daily.base, daily.basevalue.h2w,type="l",col=4)

ids2 <- a$coorx > median(a$coorx)
    
new.marks.h2e <-(a$days - as.integer(a$days))  [!ids & ids2]
temp <- hist.weighted (new.marks.h2e, wghs.daily[!ids& ids2], breaks= daily.base)
daily.basevalue.h2e <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2e <- daily.basevalue.h2e/ mean(daily.basevalue.h2e)
points(daily.base, daily.basevalue.h2e,type="l",col=5)

ids2 <- a$coory < median(a$coory)
    
new.marks.h2s <-(a$days - as.integer(a$days))  [!ids & ids2]
temp <- hist.weighted (new.marks.h2s, wghs.daily[!ids& ids2], breaks= daily.base)
daily.basevalue.h2s <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2s <- daily.basevalue.h2s/ mean(daily.basevalue.h2s)
points(daily.base, daily.basevalue.h2s,type="l",col=6)

ids2 <- a$coory > median(a$coory)
    
new.marks.h2n <-(a$days - as.integer(a$days))  [(!ids) & ids2]
temp <- hist.weighted (new.marks.h2n, wghs.daily[(!ids) & ids2], breaks= daily.base)
daily.basevalue.h2n <- ker.smooth.fft(temp$mids, temp$density, 0.03)
daily.basevalue.h2n <- daily.basevalue.h2n/ mean(daily.basevalue.h2n)
points(daily.base, daily.basevalue.h2n,type="l",col=gray(.4))

    
 legend(0.1, 1.8, lwd=2, col=c(1:6,gray(.4)),
        legend=c("All", "Downtown", "Suburb", "Suburb W",
                 "Suburb E", "Suburb S","Suburb N"),cex=0.7
        )

plot(daily.base, daily.basevalue/daily.basevalue,type="l", ylim=c(0,2), xlab="Time (days)",   ylab="Relative Level")
points(daily.base, daily.basevalue.h1/daily.basevalue,type="l",col=2)   
points(daily.base, daily.basevalue.h2/daily.basevalue,type="l",col=3)   
points(daily.base, daily.basevalue.h2w/daily.basevalue,type="l",col=4)   
points(daily.base, daily.basevalue.h2e/daily.basevalue,type="l",col=5)   
points(daily.base, daily.basevalue.h2s/daily.basevalue,type="l",col=6)   
points(daily.base, daily.basevalue.h2n/daily.basevalue,type="l",col=gray(.4))
   
 legend(0.1, 2.0, lwd=2, col=c(1:6,gray(.4)),
        legend=c("All", "Downtown", "Suburb", "Suburb W",
                 "Suburb E", "Suburb S","Suburb N"),cex=0.7
        )
  
}

if(1==1){
#x11()
###### Detect weekly periodicity difference in the first year and second years##########
ids <- a$days<=365

new.marks <- (a$days - as.integer(a$days/7)*7) 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


temp <- hist.weighted(new.marks[ids], (weights*wghs.weekly)[ids],
                      breaks=weekly.base)

weekly.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h1 <- weekly.basevalue.h1/mean(weekly.basevalue.h1)

temp <- hist.weighted(new.marks[!ids], (weights*wghs.weekly)[!ids],
                      breaks=weekly.base)

weekly.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h2 <- weekly.basevalue.h2/mean(weekly.basevalue.h2)


plot(weekly.base, weekly.basevalue, type="l",lwd=2,ylim=c(0.8,1.3), xlab="Time (days)",
     ylab="Relative Level" )
points(weekly.base, weekly.basevalue.h1, type="l",lwd=2,col=2)
points(weekly.base, weekly.basevalue.h2, type="l",lwd=2,col=3)

 legend(2, 1.2, lwd=2, col=1:3,
        legend=c("All", "2012", "2013"),cex=0.8
        )


###### Detect weekly periodicity difference in  Space: West and East##########


ids <- a$coorx<=mean(a$coorx)

new.marks <- (a$days - as.integer(a$days/7)*7) 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


temp <- hist.weighted(new.marks[ids], (weights*wghs.weekly)[ids],
                      breaks=weekly.base)

weekly.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h1 <- weekly.basevalue.h1/mean(weekly.basevalue.h1)

temp <- hist.weighted(new.marks[!ids], (weights*wghs.weekly)[!ids],
                      breaks=weekly.base)

weekly.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h2 <- weekly.basevalue.h2/mean(weekly.basevalue.h2)


plot(weekly.base, weekly.basevalue, type="l",lwd=2,ylim=c(0.8,1.3), xlab="Time (days)",
     ylab="Relative Level")
points(weekly.base, weekly.basevalue.h1, type="l",lwd=2,col=2)
points(weekly.base, weekly.basevalue.h2, type="l",lwd=2,col=3)

 legend(2, 1.2, lwd=2, col=1:3,
        legend=c("All", "West", "East"),cex=0.8
        )


###### Detect weekly periodicity difference in  Space: North and South##########


ids <- a$coory<=mean(a$coory)

new.marks <- (a$days - as.integer(a$days/7)*7) 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


temp <- hist.weighted(new.marks[ids], (weights*wghs.weekly)[ids],
                      breaks=weekly.base)

weekly.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h1 <- weekly.basevalue.h1/mean(weekly.basevalue.h1)

temp <- hist.weighted(new.marks[!ids], (weights*wghs.weekly)[!ids],
                      breaks=weekly.base)

weekly.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h2 <- weekly.basevalue.h2/mean(weekly.basevalue.h2)


plot(weekly.base, weekly.basevalue, type="l",lwd=2,ylim=c(0.8,1.3), xlab="Time (days)",
     ylab="Relative Level")
points(weekly.base, weekly.basevalue.h1, type="l",lwd=2,col=2)
points(weekly.base, weekly.basevalue.h2, type="l",lwd=2,col=3)

 legend(2, 1.2, lwd=2, col=1:3,
        legend=c("All", "South", "North"),cex=0.8
        )

###### Detect weekly periodicity difference in  Space: downtown and Suburb##########


ids = inpoly( a$coorx, a$coory, poly.downtown$V1, poly.downtown$V2)>=0

new.marks <- (a$days - as.integer(a$days/7)*7) 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


temp <- hist.weighted(new.marks[ids], (weights*wghs.weekly)[ids],
                      breaks=weekly.base)

weekly.basevalue.h1 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h1 <- weekly.basevalue.h1/mean(weekly.basevalue.h1)

temp <- hist.weighted(new.marks[!ids], (weights*wghs.weekly)[!ids],
                      breaks=weekly.base)

weekly.basevalue.h2 <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue.h2 <- weekly.basevalue.h2/mean(weekly.basevalue.h2)


plot(weekly.base, weekly.basevalue, type="l",lwd=2,ylim=c(0.8,1.4), xlab="Time (days)",
     ylab="Relative Level")
points(weekly.base, weekly.basevalue.h1, type="l",lwd=2,col=2)
points(weekly.base, weekly.basevalue.h2, type="l",lwd=2,col=3)

 legend(2, 1.2, lwd=2, col=1:3,
        legend=c("All", "Downtown", "Suburb"),cex=0.7
        )


ids2 <- a$coorx < median(a$coorx)
temp <- hist.weighted(new.marks[!ids & ids2], (weights*wghs.weekly)[!ids& ids2],
                      breaks=weekly.base)
weekly.basevalue.h2w <- ker.smooth.fft(temp$mids, temp$density, 0.5)
weekly.basevalue.h2w <- weekly.basevalue.h2w/mean(weekly.basevalue.h2w)

ids2 <- a$coorx >= median(a$coorx)
temp <- hist.weighted(new.marks[!ids & ids2], (weights*wghs.weekly)[!ids& ids2],
                      breaks=weekly.base)
weekly.basevalue.h2e <- ker.smooth.fft(temp$mids, temp$density, 0.5)
weekly.basevalue.h2e <- weekly.basevalue.h2e/mean(weekly.basevalue.h2e)

ids2 <- a$coory < median(a$coory)
temp <- hist.weighted(new.marks[!ids & ids2], (weights*wghs.weekly)[!ids& ids2],
                      breaks=weekly.base)
weekly.basevalue.h2s <- ker.smooth.fft(temp$mids, temp$density, 0.5)
weekly.basevalue.h2s <- weekly.basevalue.h2s/mean(weekly.basevalue.h2s)

ids2 <- a$coory >= median(a$coory)
temp <- hist.weighted(new.marks[!ids & ids2], (weights*wghs.weekly)[!ids& ids2],
                      breaks=weekly.base)
weekly.basevalue.h2n <- ker.smooth.fft(temp$mids, temp$density, 0.5)
weekly.basevalue.h2n <- weekly.basevalue.h2n/mean(weekly.basevalue.h2n)

plot(weekly.base, weekly.basevalue, type="l",lwd=2,ylim=c(0.8,1.4), xlab="Time (days)",
     ylab="Relative Level")
points(weekly.base, weekly.basevalue.h1, type="l",lwd=2,col=2)
points(weekly.base, weekly.basevalue.h2, type="l",lwd=2,col=3)
points(weekly.base, weekly.basevalue.h2w, type="l",lwd=2,col=4)
points(weekly.base, weekly.basevalue.h2e, type="l",lwd=2,col=5)
points(weekly.base, weekly.basevalue.h2s, type="l",lwd=2,col=6)
points(weekly.base, weekly.basevalue.h2w, type="l",lwd=2,col=gray(0.4))

 legend(2, 1.3,  lwd=2, col=c(1:6,gray(.4)),
        legend=c("All", "Downtown", "Suburb", "Suburb W",
                 "Suburb E", "Suburb S","Suburb N"),cex=0.7
        )
}

if(1==1){
#################Trend: West and East ########################################

ids <- a$coorx <= median(a$coorx)

trend.basevalue.h0 <- rep(0, length(time.marks))  
trend.basevalue.h1 <- rep(0, length(time.marks))
trend.basevalue.h2 <- rep(0, length(time.marks))

for(i in 1:nrow(a)){
     trend.basevalue.h0<- ((trend.basevalue.h0 +  wghs.trend[i]
        *dnorm(a$days[i]-time.marks, 0, 50)/
         (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
   trend.basevalue.h1 <- ((trend.basevalue.h1 + ids[i]* wghs.trend[i]*dnorm(a$days[i]-time.marks, 0, 50)/(pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
    trend.basevalue.h2 <- (trend.basevalue.h2 + (1-ids[i])* wghs.trend[i]*dnorm(a$days[i]-time.marks, 0, 50)/(pnorm(TT, a$days[i], 100)-pnorm(0, a$days[i],  50)))
}
trend.basevalue.h0 <- trend.basevalue.h0/mean(trend.basevalue.h0)

trend.basevalue.h1 <- trend.basevalue.h1/mean(trend.basevalue.h1)

trend.basevalue.h2 <- trend.basevalue.h2/mean(trend.basevalue.h2)

plot(time.marks, trend.basevalue.h0, type="l", ylim=c(0, 1.4), xlab="Time (days)",
     ylab="Relative Level")

points(time.marks, trend.basevalue.h1, type="l",col=2)

points(time.marks, trend.basevalue.h2, type="l",col=3)

 legend(80, 0.8, lwd=2, col=1:3,
        legend=c("All", "West", "East"),cex=0.8
        )

#################Trend: South and North ########################################

ids <- a$coory <= median(a$coory)
trend.basevalue.h0 <- rep(0, length(time.marks))
trend.basevalue.h1 <- rep(0, length(time.marks))
trend.basevalue.h2 <- rep(0, length(time.marks))

for(i in 1:nrow(a)){
    trend.basevalue.h0<- ((trend.basevalue.h0 +  wghs.trend[i]
        *dnorm(a$days[i]-time.marks, 0, 50)/
         (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
    trend.basevalue.h1 <- ((trend.basevalue.h1 + ids[i]* wghs.trend[i]
        *dnorm(a$days[i]-time.marks, 0, 50)/
         (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
    trend.basevalue.h2 <- (trend.basevalue.h2 + (1-ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
}
trend.basevalue.h0 <- trend.basevalue.h0/mean(trend.basevalue.h0)

trend.basevalue.h1 <- trend.basevalue.h1/mean(trend.basevalue.h1)

trend.basevalue.h2 <- trend.basevalue.h2/mean(trend.basevalue.h2)

plot(time.marks, trend.basevalue.h0, type="l", ylim=c(0, 1.4), xlab="Time (days)",
     ylab="Relative Level")

points(time.marks, trend.basevalue.h1, type="l",col=2)

points(time.marks, trend.basevalue.h2, type="l",col=3)

 legend(80, 0.8, lwd=2, col=1:3,
        legend=c("All", "South", "North"),cex=0.8
        )

############### long term trend downtown and suburb ##########
  
ids = inpoly( a$coorx, a$coory, poly.downtown$V1, poly.downtown$V2)>=0

ids.w <- a$coorx < median(a$coorx)
ids.e <- a$coorx >= median(a$coorx)
ids.s <- a$coory < median(a$coory)
ids.n <- a$coory >= median(a$coory)
        
trend.basevalue.h0 <- rep(0, length(time.marks))
trend.basevalue.h1 <- rep(0, length(time.marks))
trend.basevalue.h2 <- rep(0, length(time.marks))
trend.basevalue.h2w <- rep(0, length(time.marks))
trend.basevalue.h2e <- rep(0, length(time.marks))
trend.basevalue.h2s <- rep(0, length(time.marks))
trend.basevalue.h2n <- rep(0, length(time.marks))

for(i in 1:nrow(a)){

    trend.basevalue.h0<- ((trend.basevalue.h0 +  wghs.trend[i]
        *dnorm(a$days[i]-time.marks, 0, 50)/
         (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
    trend.basevalue.h1 <- ((trend.basevalue.h1 + ids[i]* wghs.trend[i]
        *dnorm(a$days[i]-time.marks, 0, 50)/
         (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50))))
    trend.basevalue.h2 <- (trend.basevalue.h2 + (1-ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
    trend.basevalue.h2w <- (trend.basevalue.h2w + (1-ids.w[i]*ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
    trend.basevalue.h2e <- (trend.basevalue.h2e + (1-ids.e[i]*ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
    trend.basevalue.h2s <- (trend.basevalue.h2s + (1-ids.s[i]*ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
    trend.basevalue.h2n <- (trend.basevalue.h2n + (1-ids.n[i]*ids[i])* wghs.trend[i]*
                 dnorm(a$days[i]-time.marks, 0, 50)/
                 (pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i], 50)))
    
}

trend.basevalue.h0 <- trend.basevalue.h0/mean(trend.basevalue.h0)
trend.basevalue.h1 <- trend.basevalue.h1/mean(trend.basevalue.h1)
trend.basevalue.h2 <- trend.basevalue.h2/mean(trend.basevalue.h2)
trend.basevalue.h2w <- trend.basevalue.h2w/mean(trend.basevalue.h2w)
trend.basevalue.h2e <- trend.basevalue.h2e/mean(trend.basevalue.h2e)
trend.basevalue.h2s <- trend.basevalue.h2s/mean(trend.basevalue.h2s)
trend.basevalue.h2n <- trend.basevalue.h2n/mean(trend.basevalue.h2n)

plot(time.marks, trend.basevalue.h0, type="l", ylim=c(0, 1.4), xlab="Time (days)",
     ylab="Relative Level")

points(time.marks, trend.basevalue.h1, type="l",col=2)

points(time.marks, trend.basevalue.h2, type="l",col=3)

 legend(80, 0.8, lwd=2, col=1:3,
        legend=c("All", "Downtown", "Suburb"),cex=0.8
        )


plot(time.marks, trend.basevalue.h0, type="l", ylim=c(0, 1.4), xlab="Time (days)",
     ylab="Relative Level")

points(time.marks, trend.basevalue.h1, type="l",col=2)
points(time.marks, trend.basevalue.h2, type="l",col=3)
points(time.marks, trend.basevalue.h2w, type="l",col=4)
points(time.marks, trend.basevalue.h2e, type="l",col=5)
points(time.marks, trend.basevalue.h2s, type="l",col=6)
points(time.marks, trend.basevalue.h2n, type="l",col=gray(0.4))

 legend(70, 0.8, lwd=2, col=1:3,
        legend=c("All", "Downtown", "Suburb", "Suburb W",
                 "Suburb E", "Suburb S","Suburb N"),cex=0.8
        )

    
}   ##### End of if(0==1)    

if(0==1){
############### Difference of background rate during two years ###############

background.basevalue.h1 <- matrix(0, nrow=length(background.base$x), 
   ncol=length(background.base$y))
  
background.basevalue.h2 <- matrix(0, nrow=length(background.base$x), 
   ncol=length(background.base$y))
  
if(!file.exists("Background.Smoothers"))system("mkdir Background.Smoothers")
    
for(i in 1:nrow(a)){
    print(i)
    fn <- paste("Background.Smoothers/", "bgsmoother-",i,".val", sep="")
    
   if(!file.exists(fn)){
       
       bgsmoother <- dnorm(background.basex, 
                           a$coorx[i], a$bandwidth[i])* dnorm(background.basey,
                           a$coory[i], a$bandwidth[i])/a$bg.integral[i]
       save(bgsmoother,file=fn)
    } else{
        load(fn)
    }
    
   if(a$days[i] <= 365) background.basevalue.h1 <- background.basevalue.h1 + bgprobs[i]*bgsmoother
   if(a$days[i] >= 365) background.basevalue.h2 <- background.basevalue.h2 + bgprobs[i]*bgsmoother
  }


   background.basevalue.h1 <-  background.basevalue.h1/mean(background.basevalue.h1[background.marks>=0])

   background.basevalue.h2 <-  background.basevalue.h2/mean(background.basevalue.h2[background.marks>=0])

}

if(1==1){

    ######## Plot h1
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],
                  background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue.h1*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),
                                                            seq(1, length(background.base$y), 18)],
     xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
     plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
                title(main="(a)                  "); polygon(poly.downtown, border=2)
     },
     key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),
                                              expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),
                                              expression(10^{-2}))
                    )
     } 
  )

   ########## plot h2 ##############
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue.h2*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
 plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
     title(main="(a)                  "); polygon(poly.downtown, border=2)
 },
 key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2})))} 
)


    ####### plot all ######### 
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
 plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
     title(main="(a)                  "); polygon(poly.downtown, border=2)
 },
 key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2})))} 
)

     ###### plot relative difference (h2-h1)/(h2+h1)
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (((background.basevalue.h2-background.basevalue.h1)/(background.basevalue.h2+background.basevalue.h1)*2*(background.marks>=0)/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
 plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
     title(main="(a)                  "); polygon(poly.downtown, border=2)
 },
# key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2})))} 
)

     ####### plot absolute difference h2-h1 ##############
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (((background.basevalue.h2-background.basevalue)*(background.marks>=0)/707/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
 plot.axes={axis(1); axis(2);polygon(city.boundary$x, city.boundary$y);
     title(main="(a)                  "); polygon(poly.downtown, border=2)
 },
# key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2})))} 
)

    ####### plot all with police stations #########
    police.stations <- read.csv('../police_stations.csv')
     police.stations$X = as.double(police.stations$X)/1000
      police.stations$Y = as.double(police.stations$Y)/1000
    
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], xlab="X (km)",
     ylab="Y (km)" ,
     main=(expression(mu[b](x,y))),
     plot.axes={axis(1); axis(2);
          points (police.stations$X, police.stations$Y, pch=17, cex=1.2) ;
          polygon(city.boundary$x, city.boundary$y);
          title(main="(a)                  "); polygon(poly.downtown, border=2);
 },
     key.axes={axis(4,at=seq(-9, -2), label=c(expression(10^{-9}),expression(10^{-8}),expression(10^{-7}),
                                              expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),
                                              expression(10^{-3}),expression(10^{-2})))})
    
}
dev.off()
