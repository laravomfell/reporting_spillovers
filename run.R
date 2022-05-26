# This run file reproduces all the analyses for the project
#`(No) Spillovers in domestic abuse reporting'
# by Lara Vomfell and Jan Povala.

# Author: Lara Vomfell
# Date: 11/02/2021

# The original underlying the project cannot be shared. Also, we cannot identify
# the location of the police force who provided us with the data.
# As a consequence, we also cannot share some of the GIS operations but we
# try to provide as much of the analysis code as possible.

# In the real example, we create a bounding grid over the shapefile
# and then evaluate over the entire grid.
# For the synthetic data, we generate events on a circle bounded by [0,1]
# such that the bounding grid is the unit square.

# note that running `5_plot_deprivation.R` and `6_plot_detached.R` will not work
# since in the original, we match on 2011 ONS census unit strings
# which our synthesized data does not have.

# We also note that quite a few code operations are parallelized.
# This code was written on a 64-bit Windows 10 machine, if you use a different
# operating system you might have to change some steps.

# Preliminaries ---------------------------------------------------------------

# required libraries
library(tictoc)
library(spatstat)
library(purrr)
library(Matrix)
library(fields)
library(sf)
library(data.table)
library(optparse)

library(polyCub)
library(mvtnorm)
gpclibPermit()

# set up parallel
library(foreach)
library(doParallel)
no_cores <- parallel::detectCores()

no_cores <- min(no_cores, 8)


cl <- makeCluster(no_cores)
registerDoParallel(cl)

# A few plot settings
library(snakecase)
library(ggplot2)
library(scales)
library(qqplotr)

# default plot theme
theme_print <- function(...){
  theme_classic() +
    theme(strip.text = element_text(colour = 'black'),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    theme(...)
}
theme_set(theme_print())

# define a custom labelling function for the scale to turn 0.001 into 10^-3
label_custom <- function(...){
  function(x) label_parse()(gsub("10\\^00", "1", gsub("\\+", "", gsub("1e", "10^", scientific_format()(x)))))
}
# similar idea, but also works for numbers that are not multiples of 10
label_10 <- function(...){
  function(x) label_parse()(gsub("e", "~x~10^", gsub("0e\\+00", "0", scientific_format()(x))))
}

A4_SHORT <- 8.25
A4_LONG <- 11.75

# set parameter precision
tol <- 1e-5

# set seed to jitter
set.seed(1)

# run once
if (!file.exists("poly.dll") & !file.exists("poly.so")){
  stop("You need to run 'R CMD SHLIB poly.f' on the command line once to install a (fast) Fortran version of inpoly check")
}

# import utility functions for kernels
source("utils.R")

# create folders for auxiliary files
system("mkdir background_smoothers")
system("mkdir h_space_marks")
system("mkdir figures")
system("mkdir results")

###############################################################################
#                                  MAIN PART                                  #
###############################################################################
default_start_date <- "2018-01-01"
default_end_date <- "2018-12-31"

###############################################################################
#                                   PARSING                                   #
###############################################################################

option_list = list(
  make_option("--real_data", action="store_true", default=FALSE),
  make_option("--snapshot", action="store_true", default=FALSE),
  make_option("--voronoi", action="store_true", default=FALSE),
  make_option("--g_delay", action="store_true", default=FALSE),
  make_option(c("-n", "--numsmooth"), type="integer", default=10, 
              help="Number of neighbours for background smoothing, [default= %default]", metavar="integer"),
  make_option(c("-i", "--experimentid"), type="character", default="synthetic_good",
              help="Experiment name for easier identification of files and results.", metavar="character"),
  make_option("--bw_daily", type="numeric", default=1/24.0,
              help="Bandwidth [in days] for the daily component of the background."),
  make_option("--bw_weekly", type="numeric", default=8/24.0,
              help="Bandwidth [in days] for the weekly component of the background."),
  make_option("--bw_trend", type="numeric", default=20,
              help="Bandwidth [in days] for the trend component of the background."),
  make_option("--bw_g", type="numeric", default=1.0,
              help="Bandwidth [in days] for the g(t) estimation."),
  make_option("--bw_h", type="numeric", default=0.1,
              help="Bandwidth [in km] for the h(s) estimation."),
  #
  # The following option is requried for the experiments when real_data=TRUE
  #
  make_option("--shp_fname", type="character", default="shape_file.cov"),
  #
  # Options below are there for synthetic experiemtns: real_data=FALSE
  #
  make_option("--regenerate", action="store_true", default=FALSE),
  make_option("--allevents", type="integer", default=750,
              help="Number of all events, [default= %default]", metavar="integer"),
  make_option(c("-s", "--startdate"), type="character", default=default_start_date, 
              help="Start date.", metavar="character"),
  make_option(c("-e", "--enddate"), type="character", default=default_end_date, 
              help="End date.", metavar="character"),
  make_option("--follow_trig_prob", type="numeric", default=0.1,
              help="Percentage of the triggered events that generate follow up events."),
  make_option(c("--parents_proportion"), type="numeric", default=0.8,
              help="Percentage of initial events that are considered as core. The remaining part of the initial events are triggered by a subset of core events.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


is_real_data <- opt$real_data

g_init_delay_flag <- opt$g_delay
n_p <- opt$numsmooth
bw_daily <- opt$bw_daily
bw_weekly <- opt$bw_weekly
bw_trend <- opt$bw_trend
bw_g <- opt$bw_g
bw_h <- opt$bw_h
save_snapshot <- opt$snapshot
compute_voronoi <- opt$voronoi

if (!is_real_data) {
  start_date <- as.POSIXct(opt$startdate)
  end_date <- as.POSIXct(opt$enddate)
  num_all_events <- opt$allevents
  followup_trig_prob <- opt$follow_trig_prob
  parents_proportion <- opt$parents_proportion
  regen_data <- opt$regenerate
  
  experiment_id <- paste0("synth_",
                          opt$experimentid,
                          "_numevents_", num_all_events,
                          "_np_", n_p,
                          "_bw_daily_", format(bw_daily, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_weekly_", format(bw_weekly, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_trend_", format(bw_trend, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_g_", format(bw_g, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_h_", format(bw_h, nsmall=1, digits=2, decimal.mark='_'),
                          "_follow_trig_prob_", gsub('\\.', '_', followup_trig_prob),
                          "_parents_proportion_", gsub('\\.', '_', parents_proportion),
                          "_delay_g_", tolower(g_init_delay_flag))

  # here we would normally load our shapefile, but instead we use the generated circle
  # create circle with radius 5.6km which given an area similar to the area we study.
  angle_increments <- 2 * pi / 1000
  # create angles
  ang <- seq(0, 2 * pi - angle_increments, by = angle_increments)
  circ <- data.frame(x = 5.6 + 5.6 * cos(ang),
                     y = 5.6 + 5.6 * sin(ang))
  w <- owin(poly = circ)
  shp <- st_as_sf(w)
  boundary <- data.frame(st_coordinates(shp)[, c("X", "Y")])
  boundary$x <- boundary$X
  boundary$y <- boundary$Y

  # extract the boundary
  bbox <- st_bbox(shp)

  # Synthetic data generation ---------------------------------------------------
  gen_data_id <- paste0(opt$experimentid,
                        "_numevents_", num_all_events,
                        "_np_", n_p,
                        "_follow_trig_prob_", gsub('\\.', '_', followup_trig_prob),
                        "_parents_proportion_", gsub('\\.', '_', parents_proportion))
  data_id <- gen_data_id
  gen_data_fname <- paste0("da_", gen_data_id, ".csv")
  print(paste0(">>>>>> DATA FILE NAME: ", gen_data_fname))
  
  if (file.exists(gen_data_fname) & !regen_data) {
    da <- fread(file = gen_data_fname)
  } else {
    print("Regenerating synthetic data.")
    library(stpp)
    source("generate_synthetic_data.R")
  }
} else {
  real_data_shp_fname <- opt$shp_fname
  experiment_id <- paste0("real_",
                          opt$experimentid,
                          "_np_", n_p,
                          "_bw_daily_", format(bw_daily, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_weekly_", format(bw_weekly, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_trend_", format(bw_trend, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_g_", format(bw_g, nsmall=1, digits=2, decimal.mark='_'),
                          "_bw_h_", format(bw_h, nsmall=1, digits=2, decimal.mark='_'),
                          "_delay_g_", tolower(g_init_delay_flag))
  data_id <- paste0("da_police_data_", "np_", n_p)
  
  shp <- read_sf(real_data_shp_fname, crs = 27700)
  # extract the boundary
  bbox <- st_bbox(shp) / 1000

  boundary <- data.frame(st_coordinates(shp)[, c("X", "Y")])
  boundary$X <- boundary$X / 1000
  boundary$Y <- boundary$Y / 1000
  
  # reverse x and y for owin
  boundary$x <- rev(boundary$X)
  boundary$y <- rev(boundary$Y)

  # This is necessary for spatial bandwidth computation
  w <- owin(c(bbox["xmin"], bbox["xmax"]), c(bbox["ymin"], bbox["ymax"]), poly = boundary)

  preprocessed_fname <- paste0(data_id, "_preprocessed.csv")
  if (!file.exists(preprocessed_fname)) {
    source("prepare_real_data.R")
  } else {
    da <- read.csv(preprocessed_fname)
  }
}

print(paste0(">>>>>> EXPERIMENT ID: ", experiment_id))

# Running the model ------------------------------------------------------------
source("model.R")

# Plot raw data summaries ------------------------------------------------------
source("plot_data_summaries.R")

# Plot model components -------------------------------------------------------
source("plot_model_components.R")

# Plot against detached -------------------------------------------------------
# source("plot_detached.R")

# Plot model fit --------------------------------------------------------------
source("plot_model_fit.R")

# Save expected counts according to the current model for a pre-defined space-time grid
source("expected_counts_grid.R")
