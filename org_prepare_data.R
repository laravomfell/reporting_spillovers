library(spatstat)

# THIS CAN BE SKIPPED IF type1_crime_datatable already exists

# Read Crime data
crimes <-read.csv(file="salida_puntos_llamada.csv",head=TRUE,sep=",")

crimes$hour <- substr(crimes$time,1,2)
crimes$minute <- substr(crimes$time,4,5)
crimes$second <- substr(crimes$time,7,8)

crimes$days <-(julian(as.Date(crimes$date), orig=as.Date('2012-1-1'))
               +as.double(crimes$hour)/24
               +as.double(crimes$minute)/24/60
               +as.double(crimes$second)/24/60/60 )

# crime type == 1
a<-crimes[crimes[,4]==1,]

a <- a[order(a$days),]
# add a little bit of jitter

## coordinate system?
a$coorx =a$coorx /1000+runif(nrow(a), -0.0005, 0.0005)
a$coory = a$coory /1000+runif(nrow(a), -0.0005, 0.0005)

# set ? bandwidth
a$bandwidth <- rep(0.005,nrow(a))

for(i in 1:nrow(a)){
  # calculate sq distances to all points not == i
  temp <- dist.squared(a$coorx[i], a$coory[i], a$coorx[-i], a$coory[-i])
  # keep fifth smallest distance above threshold???
  # maybe something like we need at least four points nearby??
  temp2 <- sqrt(sort(temp[temp >= 0.002^2])[5])
  # replace bandwidth if too small otherwise
  if(a$bandwidth[i] <= temp2){
    a$bandwidth[i] = temp2
  }
}


# BOUNDARY of CASTELLON
# Read boundary
city.boundary <- read.csv(file="castallon_city_boundary.csv")
city.boundary <- list(x=city.boundary$X, y= city.boundary$Y)

city.boundary$x= city.boundary$x /1000
city.boundary$y= city.boundary$y /1000

city.boundary$x <- rev(city.boundary$x)
city.boundary$y <- rev(city.boundary$y)

a$bg.integral <- rep(0, nrow(a))

w= owin(c(749, 755), c(4428, 4433), poly=city.boundary)

for(i in 1:nrow(a)){
  # calculate exact integral
  a$bg.integral[i] <- polyCub.exact.Gauss(w,
                                          mean=c(a$coorx[i], a$coory[i]),
                                          Sigma=diag(a$bandwidth[i], 2), plot=F)
  print(paste(i, 'of', nrow(a),":", a$bg.integral[i]))
  
}

write.table(a, file="type1_crime_data.table")