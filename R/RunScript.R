require(tidyverse)
require(readxl)
require(minpack.lm)
require(zoo)
require(gridExtra)
require(grid)
require(kableExtra)
require(parallel)
require(svDialogs)
# Need to detect OS or use OS independent functions
#options("scipen"=2, "digits"=3)

#source("R/InputProcessingFunctions.R")
#source("R/RTFunctions.R")
#source("R/UserFunctions.R")

#pkgload::load_code()

####### Input ############################################

processed_input <- process_input()

###################### Create plots of raw data with no adjustments ###############################
plot_list_before_baseline <- get_plots_before_baseline(processed_input)

######################### Fit the data #############################################
fits_list <- get_fits(processed_input)

###################### Create plots of response curves

rc_list <- get_rc_plots(processed_input)

####################### create plots of selected sensorgrams with fitted curves ########

plot_list <- get_fitted_plots(processed_input, fits_list)

################################ Plot all in a pdf #############################################

create_pdf(processed_input, fits_list, rc_list, plot_list)

############################# Output csv file #########################################
create_csv(processed_input, fits_list)


