library(tidyverse)
library(readxl)
library(minpack.lm)
library(zoo)
library(gridExtra)
library(grid)
library(kableExtra)
library(parallel)
library(svDialogs)
# Need to detect OS or use OS independent functions
#options("scipen"=2, "digits"=3)
#source("R/InputProcessingFunctions.R")
#source("R/RTFunctions.R")

#pkgload::load_code()

####### Input ############################################

processed_input <- process_input()

###################### Create plots of raw data with no adjustments ###############################
plot_list_before_baseline <- get_plots_before_baseline(processed_input)





fits_list <- get_fits(processed_input)


rc_list <- get_rc_plots(processed_input)

plot_list <- get_fitted_plots(processed_input, fits_list)

################################ Plot all #############################################

create_pdf(processed_input, fits_list, rc_list, plot_list)

############################# Output csv file #########################################
create_csv(processed_input, fits_list)


