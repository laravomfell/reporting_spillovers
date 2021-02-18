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
library(spatstat)
library(purrr)
library(Matrix)
library(fields)
library(sf)
library(data.table)
# needed for plotting files that will not work on synthetic data
# library(readxl)

library(polyCub)
library(mvtnorm)
gpclibPermit()

# set up parallel
library(foreach)
library(doParallel)
no_cores <- parallel::detectCores()
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# A few plot settings
library(snakecase)
library(ggplot2)
library(scales)

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

# set parameter precision
tol <- 1e-4

# run once
if (!file.exists("poly.dll")){
  stop("You need to run 'R CMD SHLIB poly.f' on the command line once to install a (fast) Fortran version of inpoly check")
}

# import utility functions for kernels
source("utils.R")

# create folders for auxiliary files
system("mkdir background_smoothers")
system("mkdir h_space_marks")
system("mkdir figures")


# Synthetic data generation ---------------------------------------------------
source("1_generate_data.R")

# Running the model -----------------------------------------------------------
source("2_model.R")

# Plot raw data summaries -----------------------------------------------------------
source("3_plot_data_summaries.R")

# Plot model components -------------------------------------------------------
source("4_plot_model_components.R")

# Plot against deprivation ----------------------------------------------------
# source("5_plot_deprivation.R")

# Plot against detached -------------------------------------------------------
# source("6_plot_detached.R")

# Plot model fit --------------------------------------------------------------
source("7_plot_model_fit.R")
