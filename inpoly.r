dyn.load('poly.so')

inpoly <-  function(x, y, px, py)
{
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

polyarea <- function(px,py)
{
  .Fortran('calcularea', as.double(px), as.double(py), as.integer(length(px)),area=as.double(0))$area
}

