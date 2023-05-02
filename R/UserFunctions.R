#' Process user input files and obtain options for fitting.
#'
#' Performs all functions selected in sample information, such as
#' automated dissociation window detection, automated concentration range, automated bulk shift detection and
#' returns a list object with the titration time series, processed sample information, all user inputs directing
#' file outputs and fitting options
#'
#'
#' @param files_directory The directory that contains the input files.
#' @param sample_sheet_path The full path to the sample information file.
#' @param data_file_path The full path to the titration data file.
#' @param output_file_path The full path where output should be stored. This directory needs to exist.
#' @param output_pdf The name of the file for the pdf output.
#' @param output_csv The name of the file for the csv output.
#' @param error_pdf The name of the file for error output.
#' @param num_cores The number of cores to use for parallel processing. The default is( the number of cores detected by `parallel::detectCores()`.
#' @param min_allowed_kd The minimum value for the dissociation constant. The default is 10^(-5).
#' @param max_iterations The maximum number of iterations for curve fitting. The default is 1000.
#' @param ptol Curve fitting parameter. If the proposed changes in parameters is smaller than this value, the optimization is considered converged. The default is 10^(-10)
#' @param ftol Curve fitting parameter. If the squared error between observed and predicted values is smaller than ftol, the optimization is considered converged. The default is 10^(-10)
#' @param min_RU_tol Minimum RU required for dissociation window detection
#' @param max_RU_tol Maximum RU required for dissociation window detection.
#'
#' @return A list.

#' @examples
#' # set up file paths for example
#'
#' files_directory <- system.file("extdata", package="htrSPRanalysis")
#' sample_sheet_path <- system.file("extdata", "sample_sheet.xlsx", package="htrSPRanalysis")
#' data_file_path <- system.file("extdata", "titration_data.xlsx", package="htrSPRanalysis")
#'
#' # process the input
#' processed_input <- process_input(files_directory = files_directory, sample_sheet_path = sample_sheet_path, data_file_path = data_file_path)
#' @export process_input

process_input <- function(files_directory = NULL,
                          sample_sheet_path = NULL,
                          data_file_path = NULL,
                          output_file_path = NULL,
                          output_pdf = NULL,
                          output_csv = NULL,
                          error_pdf = NULL,
                          num_cores = NULL,
                          min_allowed_kd = 10^(-5),
                          max_iterations = 1000,
                          ptol = 10^(-10),
                          ftol = 10^(-10),
                          min_RU_tol = 20,
                          max_RU_tol = 300){

  sample_sheet <- readxl::read_excel(sample_sheet_path)
  #identify if any flags present in the sample sheet
  flagspresent <- check_sample_sheet(sample_sheet, sample_sheet_path, files_directory)

  if (flagspresent){
     usr_msg <- "Error in sample sheet file. See Error_note_sample_sheet.csv for detailed description, then select updated file"
     stop(usr_msg)
  }

  #formatting sample_sheet data
  sample_info <- process_sample_sheet(sample_sheet)

  #import filter

  # readxl is *horribly* slow right now. Switching to openxlsx for the moment, may switch if readxl gets fixed.
  # check if this is new or old file format
  titration_data <- openxlsx::read.xlsx(data_file_path, colNames = TRUE, sheet = 1, startRow = 1, check.names = TRUE)

  ligand_and_ROI <- NULL

  if (titration_data[1,1] == "X"){
    # file may be in new format. We should skip first line
    titration_data <- openxlsx::read.xlsx(data_file_path,
                                          colNames = TRUE,
                                          sheet = 1,
                                          startRow = 2,
                                          check.names = TRUE)
    ligand_and_ROI <- openxlsx::read.xlsx(data_file_path,
                                          colNames = FALSE,
                                          sheet = 1,
                                          rows = 1,
                                          startRow = 1,
                                          check.names = TRUE)
    ligand_and_ROI <- tibble::tibble(ligand_and_ROI)

    ligand_and_ROI %>% tidyr::pivot_longer(cols = everything(), values_to = "ROI") %>%
      tidyr::separate("ROI", into = c("Name", "ROI"), sep = " \\(") %>%
      tidyr::separate("ROI", into = c("ROI", "Rest"), sep = "\\)") %>%
        dplyr::mutate(ROI = as.numeric(ROI)) %>%
        tidyr::separate(col = "Name", into = c("bracket", "Name")) %>%
        tidyr::separate(col = "Rest", into = c("Dash1","Dash2", "Analyte", "Conc String", "Conc"), sep = " ") %>%
        tidyr::separate("Conc", into = c("Conc.", "Cycle"), sep = "\\(") %>%
        dplyr::select(Name, ROI, Analyte, Conc., Cycle) %>%
        dplyr::mutate(Conc. = as.numeric(Conc.)) %>%
        dplyr::mutate(Cycle = as.numeric(Cycle)) -> ligand_and_ROI
  } else
      titration_data <- openxlsx::read.xlsx(data_file_path,
                                          colNames = TRUE,
                                          sheet = 1,
                                          startRow = 1,
                                          check.names = TRUE)

  titration_data <- tibble::tibble(titration_data)

  #identify if any flags present in the sample sheet
  flagspresent <- check_titration_data(titration_data, files_directory)

  if (flagspresent){
    usr_msg <- "Error in titration data file. See Error_note_titration_data.csv for detailed description, then select corrected file"
    stop(usr_msg)
  }


  usr_msg <- check_sample_and_data_match(sample_info, ligand_and_ROI)

  if (!is.null(usr_msg))
    stop(usr_msg)

  ################################ Set fitting options - not allowing user to change ptol or ftol at the moment ###############################

  if (is.na(min_allowed_kd) | min_allowed_kd < 10^(-7) | min_allowed_kd > 10^(-3)){
      usr_msg <- "Invalid min_allowed_kd. Please use a number in the form 1e-n, where n is between 3 and 7"
      stop(usr_msg)
  }
  if (is.na(max_iterations) | max_iterations < 500 | max_iterations > 100000){
     usr_msg <- "Invalid max_iterations. Please use a number between 500 and 100000"
     stop(usr_msg)
  }

  if (is.na(min_RU_tol) | min_RU_tol < 0 | min_RU_tol > 300){
      usr_msg <- "Invalid min_RU_tol. Please use a number between 0 and 300"
      stop(usr_msg)
  }

  if (is.na(max_RU_tol) | max_RU_tol < 50 | max_RU_tol > 500){
      usr_msg <- "Invalid max_RU_tol. Please use a number between 50 and 500"
      stop(usr_msg)
  }

  detected_num_cores <- parallel::detectCores()

  if (detected_num_cores < 1)
    detected_num_cores <- 1

  if (is.null(num_cores))
     num_cores <- detected_num_cores

  if (is.na(num_cores) | num_cores > detected_num_cores | num_cores < 1){
      usr_msg <- paste("Invalid value for num_cores. Please use a number between 1 and", detected_num_cores)
      stop(usr_msg)
  }

  ########### Set output file names #######################################
  output_file_path_default <- stringr::str_split(sample_sheet_path, "-", n = 2)[[1]][1]
  output_pdf_default <- paste0(output_file_path, date(), " - output.pdf")
  output_csv_default <- paste0(output_file_path, date()," - output.csv")
  error_pdf_default <- paste0(output_file_path, date()," - error.pdf")

  if (is.null(output_file_path))
    output_file_path <- output_file_path_default

  if (is.null(output_pdf))
    output_pdf <- output_pdf_default

  if (is.null(output_csv))
    output_csv <- output_csv_default

  if (is.null(error_pdf))
    error_pdf <- error_pdf_default

  ####### Process sample sheet #############################################
  # There are different numbers of concentrations exported for each well.
  # Here, we delete the observations from the time series and from the ligand_conc data frame
  # that are chosen by the user as not to be included.
  # For some analyses, we will restrict to the chosen concentrations.

  # keep track of ROI
  n_ROI <- dim(sample_info)[1]
  sample_info$ROI <- 1:n_ROI

  #remove wells from ligand each time series.

  selected_samples <- select_samples(sample_info, titration_data)
  #selected_samples

  expanded_sample_sheet <- sample_info
  sample_info <- selected_samples$sample_info
  keep_concentrations <- selected_samples$keep_concentrations
  Time <- selected_samples$Time
  RU <- selected_samples$RU
  all_concentrations_values <- selected_samples$all_concentrations_values
  n_time_points <- selected_samples$n_time_points

  nwells <- dim(sample_info)[1]
  # get index of first concentration (selected for analysis) for each well


  first_conc_idx_list <- purrr::map(.x = 1:nwells, .f = first_conc_indices,
                             sample_info$NumConc)

  sample_info$FirstConcIdx <- purrr::as_vector(purrr::flatten(first_conc_idx_list))

  # get baseline averages for each well (first time point)
  sample_info$BaselineAverage <- rep(0, nwells)

  baseline_info_list <- purrr::map_dfr(.x = 1:nwells, .f = get_baseline_indices,
                                sample_info,
                                Time, RU)

  sample_info$BaselineAverage <- baseline_info_list$min_baseline
  sample_info$BaselineIdx <- baseline_info_list$baseline_idx
  sample_info$BaselineNegative <- baseline_info_list$baseline_negative
  rm(baseline_info_list)

  sample_info$WellIdx <- 1:nwells

  corrected_RU <- purrr::map_dfc(.x = 1:nwells, .f = baseline_correction,
                          Time,
                          RU,
                          sample_info)


  ############################  Concentrations for fits ##############################################

  selected_concentrations <- select_concentrations(sample_info, Time, corrected_RU)

  keep_concentrations <- selected_concentrations$keep_concentrations
  sample_info <- selected_concentrations$sample_info

  incl_concentrations_values <- selected_concentrations$incl_concentrations_values
  incl_concentrations_ligand <- selected_concentrations$incl_concentrations_ligand

  error_idx_concentrations <- selected_concentrations$error_idx

  # Delete NumConc corresponding to error_idx from all FirstConcIdx after the error_idx
  # Also need to delete corresponding rows in Time and RU data frames
  # and fix all_concentrations_values

  wells <- which(!(1:nwells %in% error_idx_concentrations))
  n_fit_wells <- length(wells)

  # now the list doesn't match with sample info
  # we'll create a new sample sheet, because we need the full sheet later to display error information
  # any errors now make the first concentration index off.

  sample_info_fits <- sample_info[wells, ]

  sample_info_fits$NumInclConc <- incl_concentrations_ligand

  first_incl_conc_idx_list <- purrr::map(.x = 1:n_fit_wells, .f = first_conc_indices,
                                         sample_info_fits$NumInclConc)

  sample_info_fits$FirstInclConcIdx <- purrr::as_vector(purrr::flatten(first_incl_conc_idx_list))


  bulkshift <- purrr::map(.x = 1:n_fit_wells, .f = get_auto_bulkshift,
                              sample_info_fits,
                              Time[, keep_concentrations],
                              corrected_RU[, keep_concentrations])

  sample_info_fits$Bulkshift <- as.vector(bulkshift)
  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  parallel::clusterEvalQ(cl, library("htrSPRanalysis"))

  end_dissoc_list <- parallel::parLapply(cl, X = 1:n_fit_wells, fun = purrr::safely(find_dissociation_window),
                              sample_info_fits,
                              Time[, keep_concentrations],
                              corrected_RU[, keep_concentrations],
                              incl_concentrations_values,
                              max_RU_tol,
                              min_RU_tol)
  parallel::stopCluster(cl)

  sample_info_fits$DissocEnd <- rep(NA, n_fit_wells)

  sample_info_fits$DissocEnd <- purrr::map_dbl(.x = end_dissoc_list,
                                        .f = function(x){ifelse(is.null(x$error) & !is.null(x$result), x$result, NA)})

  list(expanded_sample_sheet = expanded_sample_sheet,
       sample_info = sample_info,
       sample_info_fits = sample_info_fits,
       Time = Time, RU = RU, corrected_RU = corrected_RU,
       keep_concentrations = keep_concentrations,
       all_concentrations_values = all_concentrations_values,
       incl_concentrations_values = incl_concentrations_values,
       n_time_points = n_time_points, max_RU_tol = max_RU_tol, min_RU_tol = min_RU_tol, nwells = nwells,
       n_fit_wells = n_fit_wells, num_cores = num_cores, min_allowed_kd = min_allowed_kd,
       max_iterations = max_iterations, ptol = ptol, ftol = ftol, output_pdf = output_pdf,
       output_csv = output_csv, error_pdf = error_pdf, error_idx_concentrations = error_idx_concentrations)
}

#' Plot all raw data that has been selected to be processed (via the `Incl.` column in the sample information). No
#' adjustments are made to the data.
#'
#' @param processed_input The list file that is output from the `process_input` function.
#' @return A list of all plots that have been selected via the `Incl.` column in sample information
#' @export get_plots_before_baseline



get_plots_before_baseline <- function(processed_input){
  sample_info <- processed_input$sample_info
  Time <- processed_input$Time
  RU <- processed_input$RU
  incl_concentrations_values <- processed_input$incl_concentrations_values
  all_concentrations_values <- processed_input$all_concentrations_values
  n_time_points <- processed_input$n_time_points
  num_cores <- processed_input$num_cores
  nwells <- processed_input$nwells

  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  parallel::clusterEvalQ(cl, library("htrSPRanalysis"))

  plot_result <- parallel::parLapply(cl, X = 1:nwells, fun = plot_sensorgrams, sample_info,
           Time, RU,
           incl_concentrations_values,
           all_concentrations_values,
           n_time_points, all_concentrations = TRUE)
  parallel::stopCluster(cl)
  plot_result
}

#' Get fits of all selected sensorgrams as indicated in the sample information.
#' @param processed_input The processed_input object returned by the function `process_input`.
#' @return A list of all fits. The fits are performed using the `safely` function, so that the list has a `$result` entry
#' and a `$error` entry for each item. If `$error` is `NULL`, the sensorgram was fit succesfully.
#' @export get_fits

get_fits <- function(processed_input){

  n_fit_wells <- processed_input$n_fit_wells
  n_time_points <- processed_input$n_time_points
  sample_info_fits <- processed_input$sample_info_fits
  num_cores <- processed_input$num_cores
  Time <- processed_input$Time
  RU <- processed_input$corrected_RU
  keep_concentrations <- processed_input$keep_concentrations
  incl_concentrations_values <- processed_input$incl_concentrations_values
  min_allowed_kd <- processed_input$min_allowed_kd
  max_iterations <- processed_input$max_iterations
  ptol <- processed_input$ptol
  ftol <- processed_input$ftol


  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  parallel::clusterEvalQ(cl, library("htrSPRanalysis"))

  fit_result <- parallel::parLapply(cl, X = 1:n_fit_wells, fun = purrr::safely(fit_association_dissociation), sample_info_fits,
           Time[, keep_concentrations],
           RU[, keep_concentrations],
           incl_concentrations_values,
           n_time_points,
           min_allowed_kd,
           max_iterations,
           ptol,
           ftol)

  parallel::stopCluster(cl)
  fit_result


}


#' Plot fitted sensorgras and raw data.
#' @param processed_input processed_input as returned by `process_input`
#' @param fits_list List of fits as returned by `get_fits`
#' @export get_fitted_plots

get_fitted_plots <- function(processed_input, fits_list){
  n_fit_wells <- processed_input$n_fit_wells
  n_time_points <- processed_input$n_time_points
  num_cores <- processed_input$num_cores
  sample_info_fits <- processed_input$sample_info_fits
  keep_concentrations <- processed_input$keep_concentrations
  incl_concentrations_values <- processed_input$incl_concentrations_values
  RU <- processed_input$corrected_RU
  Time <- processed_input$Time


  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  parallel::clusterEvalQ(cl, library("htrSPRanalysis"))

  plot_result <- parallel::parLapply(cl, X = 1:n_fit_wells, fun = plot_sensorgrams_with_fits,
           sample_info_fits, fits_list,
           Time[, keep_concentrations], RU[, keep_concentrations],
           incl_concentrations_values,
           n_time_points)

  parallel::stopCluster(cl)
  plot_result
}

#' Plot response curve. Average RU versus log10 of concentration. Color coded for concentrations selected for fitting.
#' @param processed_input Processed input object as returned from `process_input` function.
#' @export get_rc_plots

get_rc_plots <- function(processed_input){

  n_fit_wells <- processed_input$n_fit_wells
  n_time_points <- processed_input$n_time_points
  num_cores <- processed_input$num_cores
  sample_info_fits <- processed_input$sample_info_fits
  sample_info <- processed_input$sample_info
  incl_concentrations_values <- processed_input$incl_concentrations_values
  all_concentrations_values <- processed_input$all_concentrations_values
  RU <- processed_input$corrected_RU
  Time <- processed_input$Time

# Need to fix first concentration index in sample_info_fits. If we've skipped any because of errors,
# the first index is off.
  cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  parallel::clusterEvalQ(cl, library("htrSPRanalysis"))

  rc_result <- parallel::parLapply(cl, X = 1:n_fit_wells, fun = get_response_curve, sample_info_fits,
           Time, RU,
           all_concentrations_values,
           incl_concentrations_values, n_time_points)
  parallel::stopCluster(cl)
  rc_result

}

#' Create pdf file with sensorgrams with fitted curves, residuals, table of fit parameters, and response curves.
#'
#' @param processed_input#' Processed_input as returned by `process_input`
#' @param fits_list List of fits as returned by `get_fits`
#' @param rc_list List of response curves as returned by `get_rc_plots`
#' @param plot_list List of plots as returned by `get_fitted_plots`
#' @return `NULL` A pdf file is created using the path name supplied to `process_input`


#' @export create_pdf
create_pdf <- function(processed_input, fits_list, rc_list, plot_list, ...){
  nwells <- processed_input$nwells
  n_fit_wells <- processed_input$n_fit_wells
  sample_info <- processed_input$sample_info
  sample_info_fits <- processed_input$sample_info_fits
  output_pdf <- processed_input$output_pdf
  error_pdf <- processed_input$error_pdf
  error_idx_concentrations <- processed_input$error_idx_concentrations

  num_cores <- processed_input$num_cores
  #nwells <- processed_input$nwells

  # cl <- parallel::makeCluster(getOption("cl.cores", num_cores))
  # parallel::clusterEvalQ(cl, library("htrSPRanalysis"))
  #
  # pages_list <- parallel::parLapply(cl, X = 1:n_fit_wells, fun = purrr::safely(combine_output),
  #                                    fits_list,
  #                                    plot_list,
  #                                    sample_info_fits,
  #                                    rc_list,
  #                                    sample_info)
  # parallel::stopCluster(cl)
  #

  pages_list <-lapply(1:n_fit_wells, purrr::safely(combine_output),
                      fits_list, plot_list, rc_list, sample_info_fits)

  pdf(file = output_pdf, ...)
  for (well_idx in 1:n_fit_wells){
    if (!is.null(pages_list[[well_idx]]$result)){
      gridExtra::grid.arrange(pages_list[[well_idx]]$result)

    } else {
      error_msg <- paste("An error occurred when producing final output", well_idx,
                         "\n\n The sample info is: \n")
      sample_info_fits[well_idx,] %>% dplyr::select(Block, Row, Column, Ligand, Analyte) -> sample_info_error
      gridExtra::grid.arrange(grid::textGrob(error_msg), gridExtra::tableGrob(sample_info_error))
    }
  }
  dev.off()

  pdf(file = error_pdf)
  for (well_idx in 1:nwells){
    if (well_idx %in% error_idx_concentrations){
      error_msg <- paste("Too few concentrations selected/found for well", well_idx,
                         "\n\n This could be addressed by explicitly choosing concentrations to analyze. \n The sample info is: \n")
      sample_info[well_idx,] %>% dplyr::select(Block, Row, Column, Ligand, Analyte) -> sample_info_error
      gridExtra::grid.arrange(grid::textGrob(error_msg), gridExtra::tableGrob(sample_info_error))
    }
  }
  dev.off()

}

#' Create csv file with all fit parameters.
#'
#' @param processed_input#' Processed_input as returned by `process_input`
#' @param fits_list List of fits as returned by `get_fits`
#' @return `NULL` A csv file is created using the path name supplied to `process_input`


#' @export create_csv
create_csv <- function(processed_input, fits_list){

  sample_info_fits <- processed_input$sample_info_fits
  n_fit_wells <- processed_input$n_fit_wells
  output_csv <- processed_input$output_csv

  csv_data <- purrr::map_dfr(.x = 1:n_fit_wells, .f = get_csv, fits_list, sample_info_fits)

  csv_data$Ligand <- sample_info_fits$Ligand
  csv_data$Analyte <- sample_info_fits$Analyte
  csv_data$Block <- sample_info_fits$Block
  csv_data$Row <- sample_info_fits$Row
  csv_data$Column <- sample_info_fits$Column

  csv_data %>% dplyr::relocate(Ligand, Analyte, Block, Row, Column) -> csv_data

  readr::write_csv(csv_data, file = output_csv)

  csv_data
}
