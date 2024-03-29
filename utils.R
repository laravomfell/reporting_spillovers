# First step is to load function \texttt{inpoly}, which find whether each of an array of points is inside or outside of a given polygon.
if(.Platform$OS.type == "unix") {
  dyn.load('poly.so')
} else {
  dyn.load('poly')
}

inpoly <-  function(x, y, px, py){
  .Fortran("polyse",
           as.double(px), 
           as.double(py),
           as.integer(length(px)),
           as.double(x),
           as.double(y), 
           as.integer(length(x)),
           flag=integer(length(x))
  )$flag
}

polyarea <- function(px,py){
  .Fortran('calcularea', as.double(px), as.double(py), as.integer(length(px)),area=as.double(0))$area
}

# this is Simpson's integral approximation where you approx
# int_a^b f(x)dx by splitting [a,b] into m chunks and then calcul
# (b-a)/2m * (f(a) + 4f(x1) + 2f(x2) + ..+ 4(f(x_m-1) + f(b))
# This function already takes f(a), f(x1) etc and then just 
# calculates the final result
# dx gives the chunk size
simpson <- function(value, dx = 1){
  # how many chunks were given
  m <- length(value)
  
  # set up the chunk multipliers
  bx <- c(1, rep(c(4, 2), (m - 3)/2), 4, 1)
  
  sum(value * bx) * dx/3.0
}

# same procedure, just now with double integral
# as before, value is already the function evaluated (eg f(x1))
# at the chunks; this function just puts everything together
simpson.2D <- function(value, dx=1, dy=1){
    m <- nrow(value)
    n <- ncol(value)
    
    bx <- c(1, rep(c(4, 2), (m - 3)/2), 4, 1)
    by <- c(1, rep(c(4, 2), (n - 3)/2), 4, 1)
    
    sum((bx %*% value %*% by) * dx * dy/9.0)   
}

ker.smooth.conv <- function(x, y, bandwidth){
  yy <- c(y,rev(y))

  n <- length(x)
  dx <- x[2]-x[1]
  xx <- ((1:n)*dx-dx/2) 
  
  ker.val <- dnorm(xx, 0, sd=bandwidth)
  ker.val <- c(ker.val, rev(ker.val))

  smoothed <- convolve(yy, ker.val)*dx
  
  (smoothed[1:(n+1)])
}

# this function takes in a vector x and a break vector breaks
# and then essentially just calculates a histogram
# there are two tricks: they apply weights to midpoints (which they add just in case we don't observe them)
# prior to all this they also define 'tempfun' which smoothly interpolates
# between x and the breaks
hist.weighted <- function(x, weight, breaks)   {
       try(if(min(x) < min(breaks) | max(x) > max(breaks)) stop("'breaks' do not span range of 'x'") )
       nb <- length(breaks)
       
       # smooth interpolation between breaks and the length of breaks (0:200 eg)
       tempfun <- stepfun(breaks, 0:nb)
       # calculate midpoints between breaks (breaks = 0.005, ...s etc)
       mids <- (breaks[2:nb]+breaks[1:(nb-1)])*.5
       # add those midpoints to x (for smoothness?)
       x <- c(x, mids)
       # assign zero weight to those points
       weight <- c(weight, rep(0, nb-1))
       # get value of smooth interpolation at x
       z <- tempfun(x)
       # sum weights at each unique z (equal to table basically)
       freq <- tapply(weight, z, sum)
       list(breaks=breaks, 
            counts=freq, 
            density=freq/sum(freq)/(breaks[2:nb]-breaks[1:(nb-1)]), 
            mids=mids)
}


hist.weighted.2D <- function(x, y, weight, x.breaks, y.breaks){
   try(if(min(x) < min(x.breaks) | max(x) > max(x.breaks) | min(y) < min(y.breaks)
          | max(y) > max(y.breaks) ) stop("'breaks' do not span range of 'x' or 'y'!") )   
   nbx <- length(x.breaks)
   nby <- length(y.breaks)
   
   mids.x  <- (x.breaks[2:nbx]+x.breaks[1:(nbx-1)])*.5
   mids.y  <- (y.breaks[2:nby]+y.breaks[1:(nby-1)])*.5
   mids.mat.x <- mids.x %o% rep(1, nby-1)
   mids.mat.y <- rep(1, nby-1) %o%  mids.y
   xx <- c(x, mids.mat.x)
   yy <- c(y, mids.mat.y)
   weight <- c(weight, rep(0, (nbx-1)*(nby-1)))
   
   tempfun.x <- stepfun(x.breaks, 0:nbx)
   tempfun.y <- stepfun(y.breaks, 0:nby)
   z <- tempfun.x(xx) + tempfun.y(yy)*nbx
   
   freq <- matrix(tapply(weight,z, sum), nrow=nbx-1, ncol=nby-1)
   
   x1.mat <- x.breaks %o% rep(1, nby)
   y1.mat <- rep(1, nby) %o% y.breaks
   
   list(x.breaks=x.breaks, y.breaks=y.breaks, x.mids=mids.x,
        y.mids = mids.y, counts=freq, 
        density=freq/sum(freq)/(x1.mat[2:nbx,2:nby]-x1.mat[1:(nbx-1),1:(nby-1)])/(y1.mat[2:nbx, 2:nby]-y1.mat[1:(nbx-1),1:(nby-1)]))                  
}


# this function calculates the convolution of two functions f and
# f is just the density coming from our weighted histogram
# g is c(kernel of [0,1], -kernel of [0,1])
# we then calculate f x g = F^-1(F(f) * F(g))
ker.smooth.fft<- function(x, z, bandwidth){

  nx <- length(z)
  # double z
  zz  <- c(z, z)

  # get distance between midpoints
  dx <- x[2]-x[1]

  # calculate gaussian density near distance of midpoints???
  ker.val.x <- dnorm((1:nx)*dx-dx/2, 0, sd=bandwidth)
  ker.val.x <- c(ker.val.x, rev(ker.val.x))
  
  # calulate the convolution
  smoothed <- fft(fft(zz) *fft(ker.val.x), inverse=T)/nx/2*dx
  # Re returns the real part of smoothed (ignoring the imaginary part)
  Re(smoothed[1:(nx+1)])
}

ker.smooth.2D.fft<- function(x, y, z, x.bandwidth,y.bandwidth){
  nx <- nrow(z)
  ny <- ncol(z)

  zz  <- z[c(1:nx, seq(nx, 1, -1)), c(1:ny, seq(ny, 1, -1)) ] 

  dx <- x[2]-x[1]
  dy <- y[2]-y[1]
    
  ker.val.x <- dnorm((1:nx)*dx-dx/2, 0, sd=x.bandwidth)
  ker.val.y <- dnorm((1:nx)*dy-dy/2, 0, sd=y.bandwidth)
  ker.val.x <- c(ker.val.x, rev(ker.val.x))
  ker.val.y <- c(ker.val.y, rev(ker.val.y))
   
  smoothed <- fft(fft(zz) *fft(ker.val.x %o% ker.val.y), inverse=T)/nx/ny/4*dx*dy
    
  Re(smoothed[1:(nx+1), 1:(ny+1)])
}
