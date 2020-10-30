
# Now we can make the loop. 

#mu <- 0.0002
#A=0.5

for(OuterLoop in 1:40){
    
print(paste('Loop', loop, "of OuterLoop", OuterLoop))
   
################ Evaluate lambda #########################
print('############ Iterative step ########## Evaluating lambda #########################')
    

 # stop()
   
##### Reconstruction ###############################
   
   lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
   bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events

 
### 2-1. smoothing daily 
    
for(loop1 in 1:1){    
print("### 2-1. smoothing daily ")
 

wghs.daily <- daily.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events


new.marks <- a$days - as.integer(a$days)


temp <- hist.weighted (new.marks, wghs.daily, breaks= daily.base)


daily.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.03)
 
daily.basevalue <- daily.basevalue/ mean(daily.basevalue)


plot(daily.base, daily.basevalue,type="l")
    
print("### 2-2. smoothing weekly")


wghs.weekly <-  weekly.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events

new.marks <- a$days - as.integer(a$days/7)*7 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


temp <- hist.weighted(new.marks, weights*wghs.weekly, breaks=weekly.base)

weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)

plot(weekly.base, weekly.basevalue, type="l")

    
print("### 2-3. smoothing  trend")

wghs.trend <- trend.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events

trend.basevalue <- rep(0, length(time.marks))

for(i in 1:nrow(a)){
    trend.basevalue <- (trend.basevalue + wghs.trend[i]*dnorm(a$days[i]-time.marks, 0, 50)/(pnorm(TT, a$days[i], 50)-pnorm(0, a$days[i],  50)))    
}

trend.basevalue <- trend.basevalue/mean(trend.basevalue)

 #####    trend.basevalue <- rep(1, length(time.marks)) # trend can be fixed to 1 by this command


plot(time.marks, trend.basevalue, type="l")

print('################## Setting up background rate #########################')
  
    
background.basevalue <- matrix(0, nrow=length(background.base$x), 
   ncol=length(background.base$y))
  
if(!file.exists("Background.Smoothers"))system("mkdir Background.Smoothers")
    
for(i in 1:nrow(a)){ 
    fn <- paste("Background.Smoothers/", "bgsmoother-",i,".val", sep="")
    
   if(!file.exists(fn)){
       
       bgsmoother <- dnorm(background.basex, 
                           a$coorx[i], a$bandwidth[i])* dnorm(background.basey,
                           a$coory[i], a$bandwidth[i])/a$bg.integral[i]
       save(bgsmoother,file=fn)
    } else{
        load(fn)
    }
   background.basevalue <- background.basevalue + bgprobs[i]*bgsmoother
    

  if( i == nrow(a)){# print(i)
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
       (log10(background.basevalue*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], main=paste(OuterLoop, loop1, i,sep="-"))
  }
}
  
 
   background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>=0])

    
    
  ####################Re-estimate mu #####################
    
    
 trend.fun <-  approxfun(time.marks, trend.basevalue, yleft=0, yright=0)
 weekly.fun <- function(x){
     approxfun(weekly.base, weekly.basevalue,             
                       yleft=0, yright=0)(x- as.integer(x/7)*7)
 }

 daily.fun <- function(x){
     approxfun(daily.base, daily.basevalue,
                        yleft=0, yright=0)(x-as.integer(x))
 }

 background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x, 
                                           y=background.base$y, z=background.basevalue),
                                           loc=cbind(x=c(x), y=c(y))))

    
 

 bgrates.at.events.no.mu <- (trend.fun(a$days)*weekly.fun(a$days)* daily.fun(a$days)
                            *background.spatial.fun(a$coorx,a$coory))
 bgrates.at.all.no.mu <- (mean(trend.fun(time.marks)*weekly.fun(time.marks)*daily.fun(time.marks))*TT*
                          mean(background.spatial.fun(background.basex,background.basey)*background.marks)*
                          (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))

 for(kk in 1:10){
  lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
  lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A

   bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events
 
    ##### Estimate mu and A ##############

  # A =  sum(1-bgprobs)/triggers.at.all.no.A
  # mu = sum(bgprobs) / bgrates.at.all.no.mu
 
  
  #  print(paste("mu=",mu, "A=",A, "at Loop", loop1))
 }
    
} ### End of for loop1 
    
#############Reconstructing exciting component #########################
for(loop2 in 1:5){
    
excite.temporal.edge.correction <- rep(0, nrow(a))

excite.spatial.edge.correction <- rep(0, nrow(a))

    
    
for(i in 1:nrow(a)){
    excite.temporal.edge.correction[i] <-sum(excite.temporal.fun(seq(0, TT-a$days[i], 0.005)+0.6e-5))*.005
       
    temp <- Vm(paste("crime1-",substr(kk+10000,2,5), ".mark", sep=""))
          
    excite.spatial.edge.correction [i] <- simpson.2D(temp*excite.spatial.basevalue, 0.002, 0.002)
               
}

    
temporal.repetance.fun <- approxfun(excite.temporal.base, temporal.repetance, 
                                 yleft=1, yright=1)


#spatial.repetance.fun <- function(x,y){
#    temp <- interp.surface(obj=list(x=excite.spatial.base.x, y=excite.spatial.base.y,
#                                    z= spatial.repetance), loc=cbind(x=c(x), y=c(y)))
#    temp[is.na(temp)]<-0
#    temp
#}
              
   excite.wghs <- (A*excite.temporal.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]])
                   *excite.spatial.fun(a$coorx[ij.mat[,1]]-a$coorx[ij.mat[,2]],
                        a$coory[ij.mat[,1]]-a$coory [ij.mat[,2]])
                   / lambda.at.events[ij.mat[,1]])
   
   excite.temporal.series <- hist.weighted(a$days[ij.mat[,1]]-a$days[ij.mat[,2]],
                                           excite.wghs
             #/(excite.spatial.edge.correction[ij.mat[,2]]
             #*temporal.repetance.fun(a$days[ij.mat[,1]]-a$days[ij.mat[,2]]))
                                          , 
                                           breaks=excite.temporal.base)

excite.temporal.basevalue<- ker.smooth.conv(excite.temporal.series$mids, excite.temporal.series$density, bandwidth=0.2)

excite.temproal.basevalue <- excite.temporal.basevalue / simpson(excite.temporal.basevalue, 0.005)
plot(excite.temporal.series$mids, excite.temporal.series$density, pch=".",cex=2,col=2, main=paste("mu=",mu,"A=",A))
    
points(excite.temporal.base, excite.temporal.basevalue,type="l", lwd=2)

 
excite.spatial.series <-  hist.weighted.2D(dis.mat[,1], dis.mat[,2], excite.wghs
          #/(excite.temporal.edge.correction[ij.mat[,2]])
                                  , 
            x.breaks= excite.spatial.base.x, y.breaks= excite.spatial.base.y)
             

#temp<- spatial.repetance.fun(excite.spatial.series$x.mids%o%rep(1, length(excite.spatial.series$y.mids)),
#                rep(1, length(excite.spatial.series$y.mids))%o%excite.spatial.series$x.mids)
                             
 
#excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,excite.spatial.series$y.mids, 
# excite.spatial.series$density/temp, x.bandwidth=0.1, y.bandwidth=0.1) * excite.spatial.mark2
    
excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,excite.spatial.series$y.mids, 
 excite.spatial.series$density, x.bandwidth=0.1, y.bandwidth=0.1)

excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue,0.002,0.002)

filled.contour(excite.spatial.base.x[seq(1,2001,20)], excite.spatial.base.y[seq(1,2001,20)],
               excite.spatial.basevalue[seq(1,2001,20),seq(1,2001,20)], main='Kernel somoothed')
    

temp <-  hist.weighted.2D(dis.mat[,1], dis.mat[,2], excite.wghs , x.breaks= seq(-2,2,0.2), 
                          y.breaks=  seq(-2,2,0.2))
filled.contour(temp$x.mids,  temp$y.mids,temp$density, main="Rough hist")
    
    
#### Re-estimate A #########
    
 excite.temporal.fun <- approxfun(seq(0, 15, 0.005)+0.6e-12, excite.temporal.basevalue, 
                         yleft=0, yright=0)

 excite.spatial.fun <- function(x,y){
       temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                y=excite.spatial.base.y, z=excite.spatial.basevalue),
                                       loc=cbind(x=c(x), y=c(y)))
               
   #   print(temp) 
      temp[is.na(temp)]<- 0
      temp
 }
 
   
 triggers.at.events.no.A <- rep(0, nrow(a))    
 triggers.at.all.no.A <- 0

 if(1==1){ # parallel
    data.fn = paste("Loop",OuterLoop,"-",loop,".RData",sep="")
    slave.fn = paste("Loop",OuterLoop,"-",loop,".r",sep="")
    save(excite.temporal.fun, excite.spatial.fun, a, background.basex, background.basey, background.marks,Xrange, Yrange, time.marks, excite.spatial.base.x, excite.spatial.base.y,excite.spatial.basevalue,TT,file=data.fn)
    write("rm(list=ls())",file=slave.fn)
    write("suppressMessages(require('fields'))", file=slave.fn,append=T)
    write("myrank=as.integer(system('echo $MYRANK', intern=T))", file=slave.fn, append=T)
    write("nprocs=as.integer(system('echo $NPROCS', intern=T))", file=slave.fn, append=T)
    write("print(paste('Starting', myrank, 'of', nprocs))",  file=slave.fn,append=T)
    write(paste("load('", data.fn, "')",sep=""),file=slave.fn, append=T)  
    write("i1 = as.integer(nrow(a)*(myrank-1)/nprocs)+1; i2  = as.integer(nrow(a)*myrank/nprocs)",
          file=slave.fn,append=T)
    write("mytriggers.at.events.no.A = rep(0, nrow(a))",  file=slave.fn,append=T)
    write("mytriggers.at.all.no.A=0",  file=slave.fn,append=T)
    write("for(i in i1:i2){ 
           if(as.integer(i/100)*100==i)print(i)
           mytriggers.at.events.no.A <- (mytriggers.at.events.no.A +  excite.temporal.fun(a$days-a$days[i])
     *excite.spatial.fun(a$coorx-a$coorx[i],a$coory -a$coory [i])) 
    
           temp <- excite.spatial.fun(background.basex-a$coorx[i], 
                                               background.basey-a$coory[i])
    
           mytriggers.at.all.no.A <- (mytriggers.at.all.no.A+ mean(excite.temporal.fun(time.marks-a$days[i]))*TT
      *mean((background.marks>0)*temp)* (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))
       }",  file=slave.fn,append=T)
     write("save(mytriggers.at.all.no.A, mytriggers.at.events.no.A, file=paste('out-',myrank,'.image',sep=''))",
           file=slave.fn,append=T)
     write("write('Done', file=paste('out-', myrank, '.ready',sep=''))", file=slave.fn,append=T)
     write("print(paste('Done.', myrank, 'of', nprocs))",  file=slave.fn,append=T)
     write("q('no')",  file=slave.fn,append=T)
    system(paste("cat", slave.fn))
    nprocs = 20
    for(i in 1:nprocs){system(paste("MYRANK=",i," NPROCS=",nprocs, " R --slave --quiet --no-save <", slave.fn, "&",sep=""))}
    for(i in 1:nprocs){ 
        while(!file.exists(paste('out-',i,".ready",sep=""))) Sys.sleep(0.1)
        load(paste('out-', i, '.image', sep=""))
        triggers.at.events.no.A = triggers.at.events.no.A + mytriggers.at.events.no.A 
        triggers.at.all.no.A = triggers.at.all.no.A + mytriggers.at.all.no.A 
        
        unlink(paste('out-',i,".ready",sep="") )
        unlink(paste('out-',i,".image",sep="") )        
    }       
    unlink(data.fn)
 
 }
  
 
##### Estimate mu and A ##############
    
   NegLogLikehood <- function(x){
     mu <- x[1]^2
     A = x[2]^2
       
    lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
    lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A
      
    - sum(log(lambda.at.events)) + lambda.at.all
   }
    
   res.optim <- optim(par=sqrt(c(A, mu)), NegLogLikehood, control=list(trace=6))
   mu <- res.optim $par[1]^2
   A <- res.optim $par[2]^2
#   A =  sum(1-bgprobs)/triggers.at.all.no.A
#   mu = sum(bgprobs)/bgrates.at.all.no.mu

   print(paste("mu=",mu, "A=", "at Loop", loop2))

   lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
   lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A

   bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events
   
# print(c(triggers.at.events.no.A[1:10], triggers.at.all.no.A))
    
} #### End for of loop2

   print(paste("mu=",mu, "A=",A, "at OuterLoop", OuterLoop))
     save.image(file=paste('Loop-',OuterLoop, ".image", sep=""))
     
}

