# Define a couple of functions which (linearly) interpolate between x and y
trend_fun <- approxfun(time_marks, trend_basevalue, yleft=0, yright=0)

weekly_fun <- function(x){
  approxfun(weekly_base, weekly_basevalue,             
            yleft=0, yright=0)(x - as.integer(x/7)*7)
}

daily_fun <- function(x){
  approxfun(daily_base, daily_basevalue,
            yleft=0, yright=0)(x - as.integer(x))
}

# interpolate mu_b
background_fun <- function(x,y) (interp.surface(obj=list(x = background_base$x, 
                                                         y = background_base$y,
                                                         z = background_basevalue),
                                                loc=cbind(x=c(x), y=c(y))))

# initial guess for g(t)
g_fun <- approxfun(seq(0, max_t, time_units) + 0.6e-12,
                   g_basevalue, 
                   yleft=0, yright=0)

# interpolate spatial trigger component (initial guess for h(s))
h_fun <- function(x,y){
  temp <- fields::interp.surface(obj=list(x = h_base_x, 
                                  y = h_base_y, 
                                  z = h_basevalue),
                         loc=cbind(x=c(x), y=c(y))) 
  temp[is.na(temp)] <- 0
  temp
}

trigger_fun <- function(a, i){
  g_fun(a$days - a$days[i]) * h_fun(a$coorx - a$coorx[i], a$coory - a$coory[i])
}


trigger_int_fun <- function(a, time_marks, background_basex, background_basey, background_marks, constants, i){
  mean(g_fun(time_marks - a$days[i])) * constants * mean(h_fun(background_basex - a$coorx[i], background_basey - a$coory[i])[as.vector(background_marks > 0)])
}

trigger_at_all_fun <- function(i, constants){
  g_fun(time_marks - a$days[i]) * constants * mean(h_fun(background_basex - a$coorx[i], background_basey - a$coory[i])[as.vector(background_marks > 0)])
}
