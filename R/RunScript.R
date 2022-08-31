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

# Windows does not use fork. Figure out how to use multi-core on Windows
num_cores <- detectCores()/2

if (num_cores < 1)
  num_cores <- 1

####### Input Files #############################################
# select working directory
files_directory  <-  rstudioapi::selectDirectory(
  caption <- "Select Directory",
  label <- "Select"
)
#import file
sample_sheet_path <- rstudioapi::selectFile(caption = "Select the sample sheet file",
                                            filter = "All Files (*.xlsx)",
                                            existing = TRUE,
                                            path = files_directory)
sample_sheet <- read_excel(sample_sheet_path)

#identify if any flags present in the sample sheet
flagspresent <- check_sample_sheet(sample_sheet)

# requesting re-selection of the file is empty or in different format

usr_msg <- "Error in sample sheet file. See Error_note_sample_sheet.csv for detailed description, then select updated file"

while (flagspresent) {
  sample_sheet_path <- rstudioapi::selectFile(caption = usr_msg,
                                              filter = "All Files (*.xlsx)",
                                              existing = TRUE,
                                              path = files_directory)
  sample_sheet <- read_excel(sample_sheet_path)
  flagspresent <- check_sample_sheet(sample_sheet)
}

#formatting sample_sheet data
sample_info <- process_sample_sheet(sample_sheet)

#import filter
data_file_path <- rstudioapi::selectFile(caption = "Select the data file",
                                         filter = "All Files (*.xlsx)",
                                         existing = TRUE,
                                         path = files_directory)


# check if this is new or old file format

titration_data <- read_excel(data_file_path, col_names = TRUE, skip = 0, n_max = 1000)
if (titration_data[1,1] == "X")
  # file may be in new format. We should skip first line
  titration_data <- read_excel(data_file_path, col_names = TRUE, skip = 1, n_max = 1000)

#change to double - we shouldn't need to do this if file is read correctly
#titration_data  <- as.data.frame(lapply(titration_data,as.numeric))
#identify if any flags present in the sample sheet
flagspresent <- check_titration_data(titration_data)


# requesting re-selction of the file is empty or in different format

usr_msg <- "Error in titration data file. See Error_note_titration_data.csv for detailed description, then select corrected file"
while (flagspresent) {
  data_file_path <- rstudioapi::selectFile(caption = usr_msg,
                                           filter = "All Files (*.xlsx)",
                                           existing = TRUE,
                                           path = files_directory)
  titration_data <- read_excel(data_file_path, col_names = TRUE, skip = 1, n_max = 1000)
  flagspresent <- check_titration_data(titration_data)
}

################################ Set fitting options - not allowing user to change ptol or ftol at the moment ###############################
min_allowed_kd <- 10^(-5)
max_iterations <- 1000
ptol <- 10^(-10)
ftol <- 10^(-10)

min_RU_tol <- 20
max_RU_tol <- 300

flagspresent <- TRUE

while(flagspresent){
  min_allowed_kd <- dlg_input("Please enter the minimum allowed kd", print(10^(-5)))
  min_allowed_kd <- as.numeric(min_allowed_kd$res)
  flagspresent <- FALSE
  if (is.na(min_allowed_kd) | min_allowed_kd < 10^(-7) | min_allowed_kd > 10^(-3)){
    dlg_message("Please enter a number in the form 1e-n, where n is between 3 and 7")
    flagspresent <- TRUE
  }
}

flagspresent <- TRUE

while(flagspresent){
  max_iterations <- dlg_input("Please enter the maximum number of iterations for optimization algorithm", print(1000))
  max_iterations <- as.numeric(max_iterations$res)
  flagspresent <- FALSE
  if (is.na(max_iterations) | max_iterations < 500 | max_iterations > 100000){
    dlg_message("Please enter a number between 500 and 100000")
    flagspresent <- TRUE
  }
}

flagspresent <- TRUE

while(flagspresent){
  min_RU_tol <- dlg_input("Please enter the minimum RU value for choosing dissociation window", print(20))
  min_RU_tol <- as.numeric(min_RU_tol$res)
  flagspresent <- FALSE
  if (is.na(min_RU_tol) | min_RU_tol < 0 | min_RU_tol > 300){
    dlg_message("Please enter a number between 0 and 300")
    flagspresent <- TRUE
  }
}

flagspresent <- TRUE

while(flagspresent){
  max_RU_tol <- dlg_input("Please enter the maximum RU value for choosing the dissociation window", print(400))
  max_RU_tol <- as.numeric(max_RU_tol$res)
  flagspresent <- FALSE
  if (is.na(max_RU_tol) | max_RU_tol < 50 | max_RU_tol > 500){
    dlg_message("Please enter a number between 50 and 500")
    flagspresent <- TRUE
  }
}


########### Set output file names #######################################
output_file_path <- str_split(sample_sheet_path, "-", n = 2)[[1]][1]
output_pdf <- paste0(output_file_path, date(), " - output.pdf")
output_csv <- paste0(output_file_path, date()," - output.csv")
error_pdf <- paste0(output_file_path, date()," - error.pdf")



####### Process sample sheet #############################################
# There are different numbers of concentrations exported for each well.
# Here, we delete the observations from the time series and from the ligand_conc data frame
# that are chosen by the user as not to be included.
# For some analyses, we will restrict to the chosen concentrations.

#remove wells from ligand each time series.

selected_samples <- select_samples(sample_info, titration_data)

sample_info <- selected_samples$sample_info
keep_concentrations <- selected_samples$keep_concentrations
Time <- selected_samples$Time
RU <- selected_samples$RU
all_concentrations_values <- selected_samples$all_concentrations_values
n_time_points <- selected_samples$n_time_points

rm(selected_samples)

####### Baseline correction #############################################

nwells <- dim(sample_info)[1]

# get index of first concentration (selected for analysis) for each well


first_conc_idx_list <- map(.x = 1:nwells, .f = first_conc_indices,
                           sample_info$NumConc)

sample_info$FirstConcIdx <- as_vector(flatten(first_conc_idx_list))

# get baseline averages for each well (first time point)
sample_info$BaselineAverage <- rep(0, nwells)

baseline_info_list <- map_dfr(.x = 1:nwells, .f = get_baseline_indices,
                              sample_info,
                              Time, RU)

sample_info$BaselineAverage <- baseline_info_list$min_baseline
sample_info$BaselineIdx <- baseline_info_list$baseline_idx
sample_info$BaselineNegative <- baseline_info_list$baseline_negative
rm(baseline_info_list)

sample_info$WellIdx <- 1:nwells

sample_info %>% arrange(Ligand) %>% select(WellIdx) -> sort_idx
sort_idx <- as_vector(sort_idx)

#currently sort output is ignored. Need to write a function to map ROI # to desired order in output. This
# is for printing only

###################### Create plots of raw data with no adjustments ###############################
plot_list_before_baseline <- mclapply(X = 1:nwells, FUN = plot_sensorgrams, sample_info,
                                 Time, RU,
                                 incl_concentrations_values,
                                 all_concentrations_values,
                                 n_time_points, all_concentrations = TRUE, mc.cores = num_cores)

############################# Base line correction ####################################################
# baseline correction - only performed on selected concentrations


RU <- map_dfc(.x = 1:nwells, .f = baseline_correction,
                                     Time,
                                     RU,
                                     sample_info)

############################ Get Concentrations for fits ##############################################

selected_concentrations <- select_concentrations(sample_info, Time, RU)

keep_concentrations <- selected_concentrations$keep_concentrations
sample_info <- selected_concentrations$sample_info

first_incl_conc_idx_list <- map(.x = 1:nwells, .f = first_conc_indices,
                                sample_info$NumInclConc)

sample_info$FirstInclConcIdx <- as_vector(flatten(first_incl_conc_idx_list))

incl_concentrations_values <- selected_concentrations$incl_concentrations_values

error_idx_concentrations <- selected_concentrations$error_idx

########################### Plot after baseline correction ############################################
# Can be combined with fits

# Uncomment to get sensorgrams after baseline

#plot_list <- mclapply(X = 1:nwells, FUN = plot_sensorgrams, sample_info,
#                 Time[, keep_concentrations],
#                 RU[, keep_concentrations],
#                 incl_concentrations_values,
#                 all_concentrations_values,
#                 n_time_points, all_concentrations = FALSE, mc.cores = num_cores)

# Haven't corrected for duplicated values in incl_conc (only the first to be used)

####### Fits #############################################

#find best dissociation window
#fit kd
#test for bulkshift

#fits_list_kd <- map(.x = 1:nwells, .f = safely(fit_dissociation), sample_info,
#                 Time[, keep_concentrations],
#                 RU[, keep_concentrations],
#                 incl_concentrations_values)

#fits_list_fix_kd <- map(.x = 1:nwells, .f = safely(fit_association_dissociation_fix_kd), sample_info,
#                 Time[, keep_concentrations],
#                 RU[, keep_concentrations],
#                 incl_concentrations_values, fits_list_kd)



wells <- which(!(1:nwells %in% error_idx_concentrations))
n_fit_wells <- length(wells)

# now the list doesn't match with sample info... maybe need to remove those from the sample sheet?
# we'll create a new sample sheet, because we need the full sheet later to display error information

sample_info_fits <- sample_info[wells, ]

end_dissoc_list <- mclapply(X = 1:n_fit_wells, FUN = safely(find_dissociation_window), mc.cores = num_cores, sample_info_fits,
                        Time[, keep_concentrations],
                        RU[, keep_concentrations],
                        incl_concentrations_values,
                        max_RU_tol,
                        min_RU_tol)

sample_info_fits$DissocEnd <- rep(NA, n_fit_wells)

#how to extract this?


sample_info_fits$DissocEnd <- map_dbl(.x = end_dissoc_list,
                    .f = function(x){ifelse(is.null(x$error) & !is.null(x$result), x$result, NA)})


fits_list <- mclapply(X = 1:n_fit_wells, FUN = safely(fit_association_dissociation), mc.cores = num_cores, sample_info_fits,
                 Time[, keep_concentrations],
                 RU[, keep_concentrations],
                 incl_concentrations_values,
                 min_allowed_kd,
                 max_iterations,
                 ptol,
                 ftol)

rc_list <- mclapply(X = 1:n_fit_wells, FUN = get_response_curve, sample_info_fits,
               Time, RU,
               all_concentrations_values,
               incl_concentrations_values, n_time_points, mc.cores = num_cores)


plot_list_out <- mclapply(X = 1:n_fit_wells, FUN = plot_sensorgrams_with_fits,
                     sample_info_fits, fits_list,
                     Time[, keep_concentrations], RU[, keep_concentrations],
                     incl_concentrations_values,
                     n_time_points, mc.cores = num_cores)


################################ Plot all #############################################
# sort_idx provides output order. replace 1:n_fit_wells with sort_idx when function is ready

pages_list <-lapply(1:n_fit_wells, safely(combine_output),
                  fits_list, plot_list_out, rc_list, sample_info_fits)

# rm(fits_list, plot_list_out, rc_list)

pdf(file = output_pdf)
for (well_idx in 1:n_fit_wells){
  if (!is.null(pages_list[[well_idx]]$result)){
     grid.arrange(pages_list[[well_idx]]$result)

  } else {
    error_msg <- paste("An error occurred when producing final output", well_idx,
                       "\n\n The sample info is: \n")
    sample_info[well_idx,] %>% select(Block, Row, Column, Ligand, Analyte) -> sample_info_error
    grid.arrange(textGrob(error_msg), tableGrob(sample_info_error))
  }
}
dev.off()


pdf(file = error_pdf)
for (well_idx in 1:nwells){
  if (well_idx %in% error_idx_concentrations){
    error_msg <- paste("Optimal Concentration could not be determined for well", well_idx,
                       "\n\n This could be addressed by explicitly choosing concentrations to analyze. \n The sample info is: \n")
    sample_info[well_idx,] %>% select(Block, Row, Column, Ligand, Analyte) -> sample_info_error
    grid.arrange(textGrob(error_msg), tableGrob(sample_info_error))
  }
}
dev.off()

######################## Output parameters to .csv ###################################

csv_data <- map_dfr(.x = 1:n_fit_wells, .f = get_csv, fits_list, sample_info_fits)

csv_data$Ligand <- sample_info_fits$Ligand
csv_data$Analyte <- sample_info_fits$Analyte
csv_data$Block <- sample_info_fits$Block
csv_data$Row <- sample_info_fits$Row
csv_data$Column <- sample_info_fits$Column

csv_data %>% relocate(Ligand, Analyte, Block, Row, Column) -> csv_data

write_csv(csv_data, file = output_csv)
