library(foreach)
library(doParallel)
library(polyCub)
library(spatstat) # must run R with option "--max-ppsize=100000"
gpclibPermit()

# functions
dist.squared <- function(x1, y1, x2, y2) {(x1-x2)^2+(y1-y2)^2}
source("utils.R")


# THIS CAN BE SKIPPED IF type1_crime_datatable already exists

# # Read Crime data
# crimes <-read.csv(file="salida_puntos_llamada.csv",head=TRUE,sep=",")
# 
# crimes$hour <- substr(crimes$time,1,2)
# crimes$minute <- substr(crimes$time,4,5)
# crimes$second <- substr(crimes$time,7,8)
# 
# crimes$days <-(julian(as.Date(crimes$date), orig=as.Date('2012-1-1')) 
#       +as.double(crimes$hour)/24
#       +as.double(crimes$minute)/24/60
#       +as.double(crimes$second)/24/60/60 )
# 
# # crime type == 1
# a<-crimes[crimes[,4]==1,]
# 
# a <- a[order(a$days),]
# # add a little bit of jitter
# 
# ## coordinate system?
# a$coorx =a$coorx /1000+runif(nrow(a), -0.0005, 0.0005)
# a$coory = a$coory /1000+runif(nrow(a), -0.0005, 0.0005)
# 
# # set ? bandwidth
# a$bandwidth <- rep(0.005,nrow(a))
# 
# for(i in 1:nrow(a)){
#   # calculate sq distances to all points not == i
#    temp <- dist.squared(a$coorx[i], a$coory[i], a$coorx[-i], a$coory[-i])
#    # keep fifth smallest distance above threshold???
#    # maybe something like we need at least four points nearby??
#    temp2 <- sqrt(sort(temp[temp>=0.002^2])[5])
#    # replace bandwidth if too small otherwise
#    if(a$bandwidth[i] <= temp2){
#      a$bandwidth[i] = temp2
#    } 
# }
# 

# BOUNDARY of CASTELLON
# Read boundary
city.boundary <- read.csv(file="castallon_city_boundary.csv")
city.boundary <- list(x=city.boundary$X, y= city.boundary$Y)

city.boundary$x= city.boundary$x /1000
city.boundary$y= city.boundary$y /1000

city.boundary$x <- rev(city.boundary$x)
city.boundary$y <- rev(city.boundary$y)
# 
# a$bg.integral <- rep(0, nrow(a))
# 
# w= owin(c(749, 755), c(4428, 4433), poly=city.boundary)
# 
# for(i in 1:nrow(a)){
#   # calculate exact integral
#   a$bg.integral[i] <- polyCub.exact.Gauss(w, 
#                                           mean=c(a$coorx[i], a$coory[i]),
#                                           Sigma=diag(a$bandwidth[i], 2), plot=F)
#   print(paste(i, 'of', nrow(a),":", a$bg.integral[i]))
# 
# }
# 
# write.table(a, file="type1_crime_data.table")

# read data
a <- read.table("type1_crime_data.table")

# time stamps used (from where??)
time.marks <- seq(0, 730,0.005)

# TT = length of the time interval
TT <- 730

# range of coordinates
Xrange <- c(749.440, 754.264)
Yrange <- c(4428.644, 4432.570)

# a contains the crime information
# city.boundary is a list of X and Y boundary coords
# par(mfrow=c(2,2),cex=1.5)
#  
# plot(a$coorx, a$coory,col= rainbow(80)[a$days/10],xlab="X", ylab="Y",cex=0.3, main="(a)")
# points(city.boundary$x, city.boundary$y, type="l", lwd=2, col=gray(0.5))
# 
# plot(a$days, a$coory, xlab="Days", col= rainbow(80)[a$days/10],ylab="Y",cex=0.3, main="(b)")
# plot(a$coorx, a$days, xlab="X", col= rainbow(80)[a$days/10],ylab="Days",cex=0.3, main="(c)")
# 
# plot(c(0,a$days,TT), c(0:(nrow(a)+1)) , type="s", lwd=2, xlab="Days", ylab="Cum. Freq.", main="(d)")

#################### Intitial step ###################

# par(mfrow=c(3,2))

### 1. smoothing daily 

print('#### Intitial step ##### ### 1. smoothing daily ')

daily.base <- seq(0, 1,0.005)

# calculate how far into the day each event happened. yes. 1/days
new.marks <- a$days - as.integer(a$days)

temp <- hist.weighted(new.marks, rep(1, nrow(a)), breaks= daily.base)

#smooth points from midpoints density with bw 0.03
daily.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.03)
 
daily.basevalue <- daily.basevalue/ mean(daily.basevalue)

plot(daily.base,  daily.basevalue, type="l")

### 2. smoothing weekly
print('#### Intitial step ##### ### 2. smoothing weekly ')

weekly.base <- seq(0, 7, 0.005)

new.marks <- a$days - as.integer(a$days/7)*7 
weights <- 1/ (as.integer(TT/7) + (a$days - as.integer(a$days/7)*7 > TT - as.integer(TT/7)*7)) 

temp <- hist.weighted(new.marks, weights, breaks=weekly.base)


weekly.basevalue <- ker.smooth.fft(temp$mids, temp$density, 0.5)

weekly.basevalue <- weekly.basevalue/mean(weekly.basevalue)


plot(weekly.base, weekly.basevalue, type="l")

### 3. smoothing all trend
print('#### Intitial step ##### ### 3. smoothing trend ')

trend.base <- seq(0, 730,0.005)

trend.basevalue <- rep(0, length(time.marks))

for(i in 1:nrow(a)){
    trend.basevalue <- (trend.basevalue + dnorm(a$days[i]-time.marks, 0, 100)
        /(pnorm(TT, a$days[i], 100)-pnorm(0, a$days[i],  100)))    
}

trend.basevalue <- trend.basevalue/mean(trend.basevalue)

plot(time.marks, trend.basevalue, type="l")

### 3.5. #### Set up initial spatial background rate
#### For convinience, this background rate is required to be integrate to 1
#### over the whole study area

background.base <- list(x=seq(Xrange[1],Xrange[2], 0.002), 
                        y=seq(Yrange[1],Yrange[2], 0.002))

background.basevalue <- matrix(0, nrow=length(background.base$x), 
   ncol=length(background.base$y))

background.basex <- background.base$x %o% rep(1, length(background.base$y))
background.basey <- rep(1, length(background.base$x))%o% background.base$y

# for each location, check if it's outside the boundary (-1), on the boundary (0) or inside (1)
background.marks <- matrix(inpoly( background.basex ,   
            background.basey, 
            city.boundary$x, city.boundary$y), ncol=length(background.base$y))
  
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
   background.basevalue <- background.basevalue + bgsmoother
    
  if(i == nrow(a)){# print(i)
   filled.contour(background.base$x[seq(1, length(background.base$x), 18)],background.base$y[seq(1, length(background.base$y), 18)],
     (log10(background.basevalue*(background.marks>=0)/730/
               polyarea(city.boundary$x, city.boundary$y)))[seq(1, length(background.base$x), 18),seq(1, length(background.base$y), 18)], main=i)
  }
}
 
### Standardize the background so its average is 1
   background.basevalue <-  background.basevalue/mean(background.basevalue[background.marks>=0])

# The background rate is then $\mu(t, x, y) = \mu_0\, \mu_{\mathrm{trend}}(t)\mu_{\mathrm{spatial}}(x,y)\,\mu_{\mathrm{daily}}(t)\,\mu_{\mathrm{weekly}}(t) $. 
#(4) temporal and spatial components of triggering function.

print('#### Intitial step ##### ### 4.  Setting up triggering ')


 excite.temporal.basevalue <- (0.05 + seq(0,15,0.005)/10)^(-1.03)
 
 excite.temporal.basevalue <- excite.temporal.basevalue/simpson(excite.temporal.basevalue,0.005)
 
 excite.spatial.base.x <- seq(-2, 2, 0.002)
 excite.spatial.base.y<- seq(-2, 2, 0.002)
 excite.spatial.basex= excite.spatial.base.x%o% rep(1, length(excite.spatial.base.y))
 excite.spatial.basey= t(excite.spatial.base.y%o% rep(1, length(excite.spatial.base.x)))
 
 excite.spatial.basevalue = matrix(1/(abs(excite.spatial.base.x %o% rep(1, length(excite.spatial.base.y)))^2
     + abs(rep(1, length(excite.spatial.base.x))%o% excite.spatial.base.y)^2+1), ncol=2001, nrow=2001)

 excite.spatial.basevalue <- excite.spatial.basevalue/simpson.2D(excite.spatial.basevalue,0.002,0.002)

 plot(seq(0, 15, 0.005), excite.temporal.basevalue,type="l")

# now we can evaluate the lambda function 
################ Evaluated lambda #########################
print('############ Intitial step ########## Evaluating lambda #########################')

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

 excite.temporal.fun <- approxfun(seq(0, 15, 0.005)+0.6e-12, excite.temporal.basevalue, 
                         yleft=0, yright=0)

 excite.spatial.fun <- function(x,y){
      temp <- interp.surface(obj=list(x=excite.spatial.base.x, 
                                y=excite.spatial.base.y, z=excite.spatial.basevalue),
                                       loc=cbind(x=c(x), y=c(y))) 
      temp[is.na(temp)]<- 0
      temp
 }

 bgrates.at.events.no.mu <- (trend.fun(a$days)*weekly.fun(a$days)* daily.fun(a$days)
                            *background.spatial.fun(a$coorx,a$coory))
 bgrates.at.all.no.mu <- (mean(trend.fun(time.marks)*weekly.fun(time.marks)*daily.fun(time.marks))*TT*
                          mean(background.spatial.fun(background.basex,background.basey)*background.marks)*
                          (Xrange[2]-Xrange[1])*(Yrange[2]-Yrange[1]))


 triggers.at.events.no.A <- rep(0, nrow(a))    
 triggers.at.all.no.A <- 0

 if(1==1){ # parallel
     
    OuterLoop = 0
    loop =1

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
    nprocs = 16
    for(i in 1:nprocs){system(paste("MYRANK=",i," NPROCS=",nprocs, " R --slave --quiet --no-save <", slave.fn, "&",sep=""))}
    for(i in 1:nprocs){ 
        while(!file.exists(paste('out-',i,".ready",sep=""))) Sys.sleep(0.1)
        load(paste('out-', i, '.image', sep=""))
        triggers.at.events.no.A = triggers.at.events.no.A + mytriggers.at.events.no.A 
        triggers.at.all.no.A = triggers.at.all.no.A + mytriggers.at.all.no.A 
 #       plot(mytriggers.at.events.no.A)
        unlink(paste('out-',i,".ready",sep="") )
        unlink(paste('out-',i,".image",sep="") )        
    }       
    unlink(data.fn)
 }
   
 print(c(triggers.at.events.no.A[1:10], triggers.at.all.no.A))


#Update $\mu$ and $A$ for the first time. Since the likelihood function takes the form
# \begin{equation}
#   \log L = \sum_{i=1}^n \log \lambda(t_i,x_i,y_i) -\int_0^T\iint_S \lambda(t,x,y)\,\dd x\dd y\,\dd t
#  \end{equation}
# where
# \begin{equation}
#   \lambda(t)=\mu_0\mu(t,x,y) + A \sum_{i:\,t_i <t} g(t_i-t, x_i-x, y_i-y).
#   \end{equation}
# Denote $U= \int_0^T\iint_S \mu(t,x,y)\,\dd x\dd y\,\dd t$ and $G=\int_0^T\iint_S  \sum_{i:\,t_i <t} g(t-t_i, x-x_i, y-y_i)\,\dd x\dd y\,\dd t $, which are represented by \texttt{bgrates.at.all.no.mu} and \texttt{triggers.at.all.no.A}. Then the equations $\frac{\partial}{\partial \mu_0} \log L=0$ and $\frac{\partial}{\partial A} \log L=0$ 
# can be written as
# \begin{eqnarray}
#   \sum_{i=1}^n \frac{\mu(t_i,x_i,y_i)}{\lambda(t_i,x_i,y_i)} - U &=&0 \\
#   \sum_{i=1}^n \frac{\sum_{j:t_j<t_i} g(t_j-t_i,x_j-x_i,y_j-y_i)}{\lambda(t_i,x_i,y_i)} - G &=&0 \label{eq:logl.2}
#  \end{eqnarray}
# (\ref{eq:logl.2}) can be rewritten as 
# \begin{eqnarray}
#  \frac{1}{A}\sum_{i=1}^n \left[1- \frac{\mu_0 \mu(t_i,x_i,y_i)}{\lambda(t_i,x_i,y_i)} \right]- G &=&0 
#  \end{eqnarray}
# The above equations can be solve by the following iteration:
# \begin{eqnarray}
#   A^{(k+1)}&=& \frac{n-\sum_{i=1}^n \varphi_i^{(k)}}{G}\\
#   \mu_0^{(k+1)}&=& \frac{n-A^{(k+1)} G}{U}
#   \end{eqnarray}
# where $\varphi_i^{(k)} =\frac{\mu_0^{(k)}\mu(t_i,x_i,y_i)}{\mu_0^{(k)}\mu(t_i,x_i,y_i) + A^{(k)} \sum_{j:\,t_j <t_i} g(t_j-t_i, x_j-x_i, y_j-y_i)}$, represented by \texttt{bgprobs} in the program. 
A = 0.5
mu = 0.7

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

print(paste("mu=",mu, "A=", A, "at Loop", 0))


lambda.at.events <- mu * bgrates.at.events.no.mu + A * triggers.at.events.no.A
lambda.at.all <-  mu* bgrates.at.all.no.mu + A * triggers.at.all.no.A

bgprobs <- mu * bgrates.at.events.no.mu / lambda.at.events
  
#save.image('crime-excite-P1.RData')

