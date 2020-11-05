
# Now we can make the loop. 

for(OuterLoop in 1:40){
   
  # Evaluate lambda
    
  # Reconstruction
     
  lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
  bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events

  # 2-1. smoothing daily
  wghs.daily <- daily.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events
  
  new.marks <- a$days - as.integer(a$days)
  temp <- hist.weighted (new.marks, wghs.daily, breaks= daily.base)
  daily.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.03)

  daily.basevalue <- daily.basevalue/ mean(daily.basevalue)


  # 2. smoothing weekly

  wghs.weekly <-  weekly.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events

  new.marks <- a$days - as.integer(a$days/7)*7 
  weights <- 1/(as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 


  temp <- hist.weighted(new.marks, weights*wghs.weekly, breaks=weekly.base)

  weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.5)

  weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)

  # 3. smoothing  trend

  wghs.trend <- trend.fun(a$days)*background.spatial.fun(a$coorx, a$coory)/lambda.at.events

  trend.basevalue <- rep(0, length(time.marks))

  for(i in 1:nrow(a)){
    # where is the 50 coming from? wasn't it 100 before?
    trend.basevalue <- (trend.basevalue + wghs.trend[i] * dnorm(a$days[i] - time.marks, 0, 50)/
                          (pnorm(TT, a$days[i], 50) - pnorm(0, a$days[i],  50)))    
  }

  trend.basevalue <- trend.basevalue/mean(trend.basevalue)

  # Setting up background rate

  background.basevalue <- matrix(0, 
                                 nrow=length(background.base$x), 
                                 ncol=length(background.base$y))

  for(i in 1:nrow(a)){ 
    fn <- paste("Background.Smoothers/", "bgsmoother-",i,".val", sep="")
    load(fn)
    background.basevalue <- background.basevalue + bgprobs[i]*bgsmoother
  }

  # again, standardize to mean 1
  background.basevalue <- background.basevalue / mean(background.basevalue[background.marks>=0])



  # Re-estimate mu
  trend.fun <- approxfun(time.marks, trend.basevalue, yleft=0, yright=0)

  weekly.fun <- function(x){
    approxfun(weekly.base, weekly.basevalue,
              yleft=0, yright=0)(x - as.integer(x/7)*7)
  }

  daily.fun <- function(x){
    approxfun(daily.base, daily.basevalue,
              yleft=0, yright=0)(x - as.integer(x))
  }

  background.spatial.fun <- function(x,y) (interp.surface(obj=list(x=background.base$x,
                                                                   y=background.base$y, 
                                                                   z=background.basevalue),
                                                          loc=cbind(x=c(x), y=c(y))))


  bgrates.at.events.no.mu <- (trend.fun(a$days) * weekly.fun(a$days) * daily.fun(a$days) * 
                                background.spatial.fun(a$coorx,a$coory))
  
  bgrates.at.all.no.mu <- (mean(trend.fun(time.marks) * weekly.fun(time.marks) * daily.fun(time.marks)) * TT *
                           mean(background.spatial.fun(background.basex,
                                                       background.basey) * background.marks) * K)

  lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
  lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A
  
  bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events


  # Reconstructing exciting component
  for(loop2 in 1:5){

    # HM THIS DOESN't ACTUALY APPEAR TP BE USED
    excite.temporal.edge.correction <- rep(0, nrow(a))
  
    excite.spatial.edge.correction <- rep(0, nrow(a))

    for(i in 1:nrow(a)){
    excite.temporal.edge.correction[i] <-sum(excite.temporal.fun(seq(0, TT-a$days[i], 0.005)+0.6e-5))*.005
     
    # hang on. this is wrong no?
    # i think this should be i
    temp <- Vm(paste("crime1-",substr(kk+10000,2,5), ".mark", sep=""))
        
    excite.spatial.edge.correction[i] <- simpson.2D(temp*excite.spatial.basevalue, 0.002, 0.002)
             
    }

  
    excite.wghs <- (A * excite.temporal.fun(a$days[ij.mat[,1]] - a$days[ij.mat[,2]]) *
                      excite.spatial.fun(a$coorx[ij.mat[,1]] - a$coorx[ij.mat[,2]],
                                         a$coory[ij.mat[,1]]- a$coory [ij.mat[,2]]) / 
                      lambda.at.events[ij.mat[,1]])

    excite.temporal.series <- hist.weighted(a$days[ij.mat[,1]] - a$days[ij.mat[,2]],
                                            excite.wghs, 
                                            breaks=excite.temporal.base)

    # why did the kernel bandwidth change??
    excite.temporal.basevalue <- ker.smooth.conv(excite.temporal.series$mids, 
                                                 excite.temporal.series$density, 
                                                 bandwidth=0.2)

    # normalize by integral
    excite.temproal.basevalue <- excite.temporal.basevalue / simpson(excite.temporal.basevalue, 0.005)


    excite.spatial.series <- hist.weighted.2D(dis.mat[,1], 
                                              dis.mat[,2], 
                                              excite.wghs,
                                              x.breaks= excite.spatial.base.x, 
                                              y.breaks= excite.spatial.base.y)
       

    excite.spatial.basevalue <- ker.smooth.2D.fft(excite.spatial.series$x.mids,
                                                  excite.spatial.series$y.mids, 
                                                  excite.spatial.series$density, 
                                                  x.bandwidth=0.1, 
                                                  y.bandwidth=0.1)
    # normalize by integral
    excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue, 0.002,0.002)



    # Re-estimate A #########

    excite.temporal.fun <- approxfun(seq(0, 15, 0.005)+0.6e-12, 
                                     excite.temporal.basevalue, 
                                     yleft=0, yright=0)

    excite.spatial.fun <- function(x,y){
      temp <- interp.surface(obj=list(x = excite.spatial.base.x, 
                                      y = excite.spatial.base.y, 
                                      z = excite.spatial.basevalue),
                             loc = cbind(x=c(x), y=c(y)))
      temp[is.na(temp)] <- 0
      temp
    }


    triggers.at.events.no.A <- rep(0, nrow(a))    
    triggers.at.all.no.A <- 0

    # this loop is very slow
    for (i in 1:nrow(a)){
      if (i %% 100 == 0) print(paste("on:", i))
      t_temp <- excite.temporal.fun(a$days - a$days[i]) * excite.spatial.fun(a$coorx - a$coorx[i], a$coory - a$coory[i])
      triggers.at.events.no.A <- triggers.at.events.no.A + t_temp
      
      temp <- excite.spatial.fun(background.basex - a$coorx[i], 
                                 background.basey - a$coory[i])
      
      triggers.at.all.no.A <- triggers.at.all.no.A + mean(excite.temporal.fun(time.marks - a$days[i])) * TT * 
        mean(temp[background.marks > 0]) * bg_weight * K
    }


    # Estimate mu and A ##############
    
    res.optim <- optim(par=sqrt(c(A, mu)), NegLogLikehood, control=list(trace=6))
    mu <- res.optim$par[1]^2
    A <- res.optim$par[2]^2
    
    lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
    lambda.at.all <-  mu * bgrates.at.all.no.mu + A * triggers.at.all.no.A
    
    bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events
    
    } #### End for of loop2

}

