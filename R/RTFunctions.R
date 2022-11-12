#' Automatically determine dissociation window
#'
#' @param well_idx The corresponding well in the extended sample sheet
#' @param sample_info A tibble. The expanded sample sheet
#' @param x_vals A tibble. Time columns from the carterra output for concentrations chosen for fitting
#' @param y_vals A tibble RU values from the carterra output for concentrations chosen for fitting
#' @param incl_concentrations_ligand A vector. The concentrations chosen for inclusion in fitting.
#' @param max_RU_tol A number. The maximum RU value to be included for determining the dissociation window
#' @param min_RU_tol A number. The minimum RU value to be included for determining the dissociation window
#' @return A number. The end dissociation time required to capture the decay in all of the selected concentrations for the given well
#'
#' find_dissociation_window(well_idx, sample_info, x_vals, y_vals, incl_concentrations_ligand, max_RU_tol, min_RU_tol)

find_dissociation_window <- function(well_idx, sample_info, x_vals, y_vals,
                                     incl_concentrations_ligand, max_RU_tol, min_RU_tol){

  if (sample_info[well_idx,]$`Automate Dissoc. Window` != "Y")
    return(NULL)

  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association

  start_time <- baseline + baseline_start + association
  dissociation <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1

  end_dissoc_list <- NULL
  df_all <- NULL
  max_idx <- 0
  for (i in start_idx:end_idx){
    Time <- x_vals[, i]
    RU <- y_vals[, i]
    df <- suppressMessages(dplyr::bind_cols(Time, RU))
    names(df) <- c("Time", "RU")
    df %>% dplyr::filter(Time > (start_time - 50)) -> df

    # Do not base on low information concentrations
    if (mean(df$RU, na.rm = TRUE) < min_RU_tol | mean(df$RU, na.rm = TRUE) > max_RU_tol)
      next
    df_zoo <- zoo::as.zoo(df)
    #  max_idx <- max_idx + 1
    tibble::as_tibble(
      zoo::rollapply(
        data = df_zoo, FUN =  function(x){
          x_df <- tibble::as_tibble(x)
          if (all(is.na(x_df$RU)| is.nan(x_df$RU)))
            return(c(NA,NA))
          else
            return(stats::coef(stats::lm(RU ~ Time, singular.ok = TRUE,
                           data = x_df)))}, by.column = FALSE,

        width = 100)) -> df_out
    names(df_out) <- c("Intercept", "Slope")
    n_vals <- dim(df_out)[1]

    #  df_out %>% dplyr::mutate(x = rep(max_idx, n_vals)) -> df_out
    df_out %>% dplyr::mutate(RollIndex = 1:n_vals) -> df_out
    max_slope <- max(abs(df_out$Slope))
    target_slope <- .01*max_slope
    window_idx <- which(abs(df_out$Slope) < target_slope)[1]

    if (is.na(window_idx)){
      # Once this happens, we are using the entire time series for all concentrations
      end_dissoc_list <- c((df$Time)[length(df$Time)], end_dissoc_list)
      next
    } else
         end_dissoc <- df$Time[window_idx]

    end_dissoc_list <- c(end_dissoc, end_dissoc_list)

  }
  # Now we have a candidate end of dissoc for each concentration. Thia should be done for selected concentrations
  # Overall window is smallest window that accomodates all the concentrations
  return(max(as.numeric(end_dissoc_list), na.rm = TRUE))
}

#' Determine the first concentration index in the Carterra output for a sample in the extended sample sheet
#'
#' @param well_idx The corresponding well in the extended sample sheet
#' @param num_conc_ligand A vector. For each ligand (well), the total number of analyte concentrations observed.
#' @return A number. The first index in the Carterra ouput that corresponds to this ligand (well)
#'
#' first_conc_indices(well_idx, num_conc_ligand)


first_conc_indices <- function(well_idx, num_conc_ligand){

  #Compute the correction to baseline. Usually for regenerative case.

  # find displacement from last ligand
  if (well_idx == 1){
    first_conc_idx <- 1
  } else
    first_conc_idx <- sum(num_conc_ligand[1:(well_idx - 1)]) + 1 #number of time series up to this well
  first_conc_idx
}

#' Determine the index of the baseline concentration
#'
#' @param well_idx The corresponding well in the extended sample sheet
#' @param sample_info A tibble. The expanded sample sheet
#' @param x_vals A tibble. Time columns from the Carterra output for concentrations chosen for fitting
#' @param y_vals A tibble RU values from the Carterra output for concentrations chosen for fitting
#' @param sample_info A tibble. This is the sample sheet after it has been extended to include a row for each well
#' @return A list.
#' @param baseline_idx A number. The index of the minimum average baseline for this well.
#' @param min_baseline A number. The average of the minimum baseline
#' @param baseline_negative A logical. TRUE if the baseline for the highest concentrate has a negative average.
#'
#' get_baseline_indices(well_idx, sample_info, x_vals, y_vals)

get_baseline_indices <- function(well_idx, sample_info, x_vals, y_vals){
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1

  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  baseline_avg_list <- NULL

  for (i in start_idx:end_idx){
      Time <- x_vals[, i]
      RU <- y_vals[, i]
      df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU))
      colnames(df) <- c("Time", "RU")
      df %>% dplyr::filter(Time > baseline_start & Time < baseline+ baseline_start)  %>% .$RU -> base_meas
      baseline_avg <- mean(base_meas, na.rm = TRUE)
      baseline_avg_list <- c(baseline_avg_list, baseline_avg)
  }
  min_baseline <- min(baseline_avg_list)
  baseline_idx <- start_idx + which(baseline_avg_list == min_baseline) - 1

  #is highest baseline negative?
  if (baseline_avg < 0)
    baseline_neg <- TRUE else
    baseline_neg <- FALSE

  list(baseline_idx = baseline_idx, min_baseline = min_baseline, baseline_negative = baseline_neg)
}

# internal function. Not documented
create_dataframe_with_conc <- function(begin_conc_idx, end_conc_idx, x_vals, y_vals,
                                 numerical_concentrations,
                                n_time_points){
    n_vals <- dim(x_vals)[2]
    names(x_vals) <- as.character(1:n_vals)
    names(y_vals) <- as.character(1:n_vals)

    Time <- x_vals[, begin_conc_idx:end_conc_idx] %>%
      tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
      dplyr::select(value)
    RU <- y_vals[, begin_conc_idx:end_conc_idx] %>%
      tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
      dplyr::select("value")

    purrr::map_dfr(.x = tibble::tibble(numerical_concentrations), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
      dplyr::arrange(numerical_concentrations) -> numerical_concentrations
    Concentrations <- numerical_concentrations
    df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))
    colnames(df) <- c("Time", "RU", "Concentration")
    df
}

#This function is under development
bulkshift_correction <- function(well_idx, x_vals, y_vals, sample_info,
                                 all_concentrations_ligand){

  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  association <- sample_info[well_idx,]$Dissociation
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- all_concentrations_ligand[well_idx]

}

#' Determine baseline correction
#'
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param x_vals A tibble. Time columns from the Carterra output for concentrations chosen for fitting.
#' @param y_vals A tibble RU values from the Carterra output for concentrations chosen for fitting.
#' @param sample_info A tibble. The expanded sample sheet.
#' @return A tibble. The baseline-adjusted RU values.
#'
#' baseline_correction <- function(well_idx, x_vals, y_vals, sample_info)

baseline_correction <- function(well_idx, x_vals, y_vals, sample_info){

  negative_baseline <- sample_info[well_idx, ]$BaselineNegative
  baseline_average <- sample_info[well_idx, ]$BaselineAverage
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`

  start_idx <- sample_info[well_idx, ]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc  #number of concentrations for this well
  end_idx <- start_idx + num_conc - 1


  #Need to test if baseline of highest conc is < 0

  if (sample_info[well_idx,]$Regen. == "N" & !negative_baseline){
    # add baseline average to all timepoints

    y_vals[, start_idx:end_idx] <-
      y_vals[, start_idx:end_idx] - baseline_average
  } else
    {
    # correct all to mean of zero
    for (i in 1:num_conc){
      # compute baseline average for each concentration
      # subtract from baseline average from all RU vals for that concentration

      Time <- x_vals[, start_idx + (i-1)]
      RU <- y_vals[, start_idx + (i-1)]

      df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU))
      colnames(df) <- c("Time", "RU")

      df %>% dplyr::filter(Time > baseline_start & Time < baseline + baseline_start) %>% .$RU -> base_meas # select for defined baseline time period
      # This command is split up because mean(.$RU) would not parse properly
      base_corr <- mean(base_meas, na.rm = TRUE)
      y_vals[, start_idx + (i-1)] <- y_vals[, start_idx + (i-1)] - base_corr
    }
  }
  y_vals[, start_idx:end_idx]

}

#' Obtain the response curve of log concentrations vs. RU averaged over 5 time points at the end of the association phase.
#'
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param sample_info A tibble. The expanded sample sheet.
#' @param x_vals A tibble. Time columns from the Carterra output for all concentrations.
#' @param y_vals A tibble RU values from the Carterra output for all concentrations.
#' @param all_concentrations_values A vector. All concentrations measured.
#' @param incl_concentrations_values A vector. Concentrations to include in fitting.
#' @param n_time_points A number. The maximum number of time points collected in this experiment. It is the number of rows in the Time and RU tibbles.
#' @return A ggplot2::ggplot object. The plot of the response curve.
#'
#' get_response_curve <- function(well_idx, sample_info, x_vals, y_vals, all_concentrations_values, incl_concentrations_values, n_time_points)

get_response_curve <- function(well_idx, sample_info, x_vals, y_vals,
                               all_concentrations_values,
                               incl_concentrations_values, n_time_points){

  start_incl_idx <- sample_info[well_idx,]$FirstInclConcIdx
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  end_idx <- start_idx + num_conc - 1

  num_incl_conc <- sample_info[well_idx,]$NumInclConc
  end_incl_idx <- start_incl_idx + num_incl_conc - 1
  incl_conc_values <- incl_concentrations_values[start_incl_idx:end_incl_idx]

  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association

  df <- create_dataframe_with_conc(start_idx, end_idx, x_vals, y_vals,
                                         all_concentrations_values[start_idx:end_idx],
                                         n_time_points)

  df %>% dplyr::filter((Time >= baseline + baseline_start + association - 10)
                & (Time <= baseline + baseline_start + association - 5)) %>%
    dplyr::group_by(Concentration) %>%
    dplyr::summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC
  df_RC %>% dplyr::mutate(Included =
                     forcats::as_factor(ifelse(Concentration %in% incl_conc_values,
                                      "Yes", "No"))) -> df_RC
  ggplot2::ggplot(df_RC, ggplot2::aes(x = Concentration,
                    y = AverageRU)) +
    ggplot2::geom_point(ggplot2::aes(color = Included)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle(ligand_desc)
}

#' Determine the best range of concentrations for fitting. The algorithm will attempt to find 5 concentrations
#' by computing a rolling slope of the response curve over
#' 5 consecutive concentrations. The concentration that leads to the maximum cumulative sum of slopes is marked as
#' the starting concentration. A concentration is removed from consideration if its slope (the window where it is the first concentration)
#' is less than 20% of the average of all slopes. If none of the concentrations are removed, the algorithm returns
#' the five concentrations. Otherwise, it removes the concentration with the smallest slope and returns 4 concentrations.
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param sample_info A tibble. The expanded sample sheet.
#' @param x_vals A tibble. Time columns from the Carterra output for all concentrations.
#' @param y_vals A tibble RU values from the Carterra output for all concentrations.
#' @param num_conc A number. The number of concentrations observed for this well.
#' @param concentrations A vector. The concentrations observed for this well.
#' @param start_idx A number. The start index for the columns in the titration data for the current well.
#' @return A vector. The concentrations selected for fitting.
#'
#' get_best_window <- function(well_idx, sample_info, x_vals, y_vals, num_conc, concentrations, start_idx)

get_best_window <- function(well_idx, sample_info, x_vals, y_vals,
                            num_conc, concentrations, start_idx, n_time_points){

  association <- sample_info[well_idx,]$Association
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  end_assoc_time <- association + baseline + baseline_start

  end_idx <- start_idx + num_conc - 1
  end_assoc_resp <- NULL

  association_end <- baseline + baseline_start + association

  df <- create_dataframe_with_conc(start_idx, end_idx, x_vals, y_vals,
                                         concentrations,
                                         n_time_points)

  #this has failed in some data sets where these observations are missing.
  df %>% dplyr::filter((Time >= baseline + baseline_start + association - 10)
                & (Time <= baseline + baseline_start + association - 5)) %>%
    dplyr::group_by(Concentration) %>%
    dplyr::summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC

  # check to see if any results for df_RC

  if (dim(df_RC)[1] == 0)
     return(NULL)

  # record the differences between consecutive responses
  sum_diff <- NULL
  for (i in 1:(num_conc - 1)){
    sum_diff <- suppressMessages(dplyr::bind_cols(sum_diff, df_RC[i+1,]$AverageRU - df_RC[i,]$AverageRU))
  }
  # add sums for each 5 cycle window.
  cum_sum <- zoo::rollapply(purrr::as_vector(purrr::flatten(sum_diff)), 4, FUN = sum)
  start_conc_idx <- which(cum_sum == max(cum_sum))

  # check start and end slopes, may be better to fit 4 instead of five

  slopes <- purrr::as_vector(sum_diff[start_conc_idx:(start_conc_idx+3)])

  remove_concentration <- ifelse(slopes < 0.2*mean(slopes), 1, 0)

  # return 5 best consecutive concentrations
  if(sum(remove_concentration) == 0)
     return(concentrations[start_conc_idx:(start_conc_idx + 4)])

  # if some slopes are less than 20% of the mean, remove one concentration
  # low end or high end, depending on which slope is smaller

  # remove_concentration has at least one '1' value
  # if both first and last are taggplot2::gged, we remove the smallest
  # if only one is taggplot2::gged, it will still be the smallest
  if (slopes[1] < slopes[4])
      return(concentrations[(start_conc_idx+1):(start_conc_idx + 4)])
  else
      return(concentrations[(start_conc_idx):(start_conc_idx + 3)])

}

#' Plot sensorgrams. This function will plot only the data. For plotting data with fitted curves, use plot_sensorgrams_with_fits
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param sample_info A tibble. The expanded sample sheet.
#' @param x_vals A tibble. Time columns from the Carterra output for all concentrations.
#' @param y_vals A tibble RU values from the Carterra output for all concentrations.
#' @param incl_conc_values A vector. The number of concentrations to be included in this plot.
#' @param all_concentrations_values A vector. The concentrations observed for this well.
#' @param n_time_points A number. The maximum number of time points observed in this experiment. It is
#' the number of rows in the titration data tibble.
#' @param all_concentrations A logical. Whether or not to include all observed concentrations in the plot.
#' The default is FALSE. This will plot only concentrations included in the fit.
#' @return A ggplot2::ggplot object. The plotted sensorgram.
#'
#' plot_sensorgrams <- function(well_idx, sample_info, x_vals, y_vals,
#' incl_conc_values, all_concentrations_values, n_time_points,
#' all_concentrations = FALSE)


plot_sensorgrams <- function(well_idx,
                             sample_info,
                             x_vals,
                             y_vals,
                             incl_conc_values,
                             all_concentrations_values,
                             n_time_points,
                             all_concentrations = FALSE){
  if (!all_concentrations){
    start_idx <- sample_info[well_idx,]$FirstInclConcIdx
    num_conc <- sample_info[well_idx,]$NumInclConc
  } else
  {
    start_idx <- sample_info[well_idx,]$FirstConcIdx
    num_conc <- sample_info[well_idx,]$NumConc
    incl_conc_values <- all_concentrations_values
  }
  end_idx <- start_idx + num_conc - 1

  incl_conc_values <- incl_conc_values[start_idx:end_idx]

  ligand_desc <- sample_info[well_idx,]$Ligand

  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)


  Time <- x_vals[, start_idx:end_idx] %>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select(value)
  RU <- y_vals[, start_idx:end_idx]%>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select("value")


  purrr::map_dfr(.x = tibble::tibble(incl_conc_values), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    dplyr::arrange(incl_conc_values) -> incl_conc_values

   Concentrations <- incl_conc_values

  df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))

  colnames(df) <- c("Time", "RU", "Concentration")

  df$Concentration <- forcats::as_factor(formatC(df$Concentration, format = "e",digits = 2))


  sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row", sample_info[well_idx,]$Row,
                     "Column", sample_info[well_idx,]$Column)

  ggplot2::ggplot(df, ggplot2::aes(x = Time, y = RU, color = Concentration)) + ggplot2::geom_point(size = 0.5) +
    ggplot2::ggtitle(ligand_desc, subtitle = sub_title)
}

#' Plot sensorgrams. This function will plot only the data. For plotting data with fitted curves, use plot_sensorgrams_with_fits
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param sample_info A tibble. The expanded sample sheet.
#' @param fits A kinetics fit object. This is returned from a 'safely' call to fit_association_dissociation.
#' @param x_vals A tibble. Time columns from the Carterra output for all concentrations.
#' @param y_vals A tibble RU values from the Carterra output for all concentrations.
#' @param incl_conc_values A vector. The number of concentrations to be included in this plot.
#' @param n_time_points A number. The maximum number of time points observed in this experiment. It is
#' the number of rows in the titration data tibble.
#' @return A ggplot2::ggplot object. The plotted sensorgram with fitted curves.
#'
#' plot_sensorgrams_with_fits <- function(well_idx, sample_info, fits, x_vals, y_vals,
#' incl_conc_values, n_time_points)

plot_sensorgrams_with_fits <- function(well_idx, sample_info, fits, x_vals, y_vals,
                             incl_conc_values, n_time_points){

  if (!is.null(fits[[well_idx]]$error))
    return(NULL)

  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1

  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association
  dissociation <- sample_info[well_idx,]$Dissociation

  assoc_start <- baseline + baseline_start

  fit_RU <- fits[[well_idx]]$result$FitOutcomes$RU

  n_vals <- dim(x_vals)[2]
  names(x_vals) <- 1:n_vals
  names(y_vals) <- 1:n_vals


  Time <- x_vals[, start_idx:end_idx] %>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select(value)
  RU <- y_vals[, start_idx:end_idx]%>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select("value")

  numerical_concentration <- incl_conc_values[start_idx:end_idx]

  purrr::map_dfr(.x = tibble::tibble(numerical_concentration), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    dplyr::arrange(numerical_concentration) -> numerical_concentration

  Concentrations <- numerical_concentration

  df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))

  colnames(df) <- c("Time", "RU", "Concentration")

  end_time <- ifelse(is.na(sample_info[well_idx, ]$DissocEnd), baseline + baseline_start + association + dissociation, sample_info[well_idx, ]$DissocEnd)

  df %>% dplyr::filter(Time > baseline + baseline_start & Time < end_time) -> df

  suppressMessages(dplyr::bind_cols(df,FittedRU = fit_RU)) -> df

  # Correct to zero time

  df$Time <- df$Time - assoc_start

  colnames(df) <- c("Time", "RU", "Concentration", "FittedRU")

  df$Concentration <- forcats::as_factor(formatC(df$Concentration, format = "e",digits = 2))

  sub_title <- paste("Block", sample_info[well_idx,]$Block, "Row",
                     sample_info[well_idx,]$Row,
                     "Column", sample_info[well_idx,]$Column)


  ggplot2::ggplot(df, ggplot2::aes(x = Time, y = RU)) + ggplot2::geom_point(size = 0.09, ggplot2::aes(color = Concentration)) +
    ggplot2::ggtitle(ligand_desc, subtitle = sub_title) +
    ggplot2::geom_line(ggplot2::aes(x = Time, y = FittedRU, group = Concentration), color = "black")
}

# This function is used internally. It is a replacement for summary.nlm that allows for a non-singular
# Hessian due to the constrained fits.

summary_fit_with_constraints <- function(fit_object){

  info <- fit_object$info
  hessian <- fit_object$hessian
  pars <- fit_object$par
  n <- length(pars)
  std_err_full <- rep(NA, n)

  if (info == 0 | info == 5) {
    df <- data.frame(Estimate = pars, "Std. Error" = std_err_full)
    colnames(df) <- c("Estimate", "Std. Error")
    return(df)
  }

  n <- nrow(hessian)
  test_zeroes <- apply(hessian, 1 , function(x) sum(x==0))
  nonsingular_rows <- which(test_zeroes != n)
  if (length(nonsingular_rows) == n & info != 5)
    return(summary(fit_object)$coefficients)
  hessian <- hessian[nonsingular_rows, nonsingular_rows]

  std_err_full <- rep(NA, n)

# get table directly. Code is pulled from summary.minpack.lm

  if (info != 5) {            # when info is 5, that means the iterations maxed out and fit is not valid
    ibb <- chol(hessian)
    ih <- chol2inv(ibb)
    p <- length(pars)
    rdf <- length(fit_object$fvec) - p
    resvar <- stats::deviance(fit_object)/rdf
    se <- sqrt(diag(ih) * resvar)
  }
  else
    se <- rep(NA, length(nonsingular_rows))

  std_err_full[nonsingular_rows] <- se

  df <- data.frame(Estimate = pars, "Std. Error" = std_err_full)
  colnames(df) <- c("Estimate", "Std. Error")
  df
}

# Fit only kd. Not implemented
fit_kd <- function(pars, df, incl_concentrations, num_conc, kd, t0 = t0){
  #pars ("Rmax" one for each concentration,"ka", "tstart" one for each concentration)

  R0 <- pars[1:num_conc]
  kd <- pars[(num_conc+1)]

  err_assoc <- NULL
  err_dissoc <- NULL

  for (i in 1:num_conc){

    df_i <- df %>% dplyr::filter(Concentration == incl_concentrations[i])
    RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration

    df_i %>% dplyr::filter(DissocIndicator == 1) -> df_dissoc

    dissoc_formula <- R0[i]*exp(-kd*(df_dissoc$Time - t0))

    err_dissoc <-   c(err_dissoc, df_dissoc$RU - dissoc_formula)

  }
  err_dissoc

}

# Internal function that is passed to nlm. It computes the objective function for the fit.
fit_as_system <- function(pars, df, incl_concentrations, num_conc, association, bulkshift, global_rmax){

  if(global_rmax){
    #pars (global "Rmax","ka", "R0" one for each concentration, kd, shift one for each concentration)

    Rmax <- pars[1]

    ka <- pars[2]
    t0 <- pars[3:(2 + num_conc)]
    kd <- pars[(3+num_conc)]

    if (bulkshift)
      shift <- pars[(4 + num_conc):(3 + 2*num_conc)]

  } else{
    #pars ("Rmax" one for each concentration,"ka", "R0" one for each concentration, kd, shift one for each concentration)

    Rmax <- pars[1:num_conc]

    ka <- pars[num_conc+1]
    t0 <- pars[(num_conc + 2):(2*num_conc + 1)]
    kd <- pars[2*num_conc + 2]

    if (bulkshift)
      shift <- pars[(2*num_conc + 3):(3*num_conc + 2)]

  }

  err_assoc <- NULL
  err_dissoc <- NULL

  for (i in 1:num_conc){

    df_i <- df %>% dplyr::filter(Concentration == incl_concentrations[i])
  #  RU <- df_i$RU

    df_i %>% dplyr::filter(AssocIndicator == 1) -> df_assoc
    df_i %>% dplyr::filter(DissocIndicator == 1) -> df_dissoc

    if (global_rmax){
      assoc_formula_first_term <-
        (Rmax * ka * incl_concentrations[i])/(ka*incl_concentrations[i] + kd)

    } else {
      assoc_formula_first_term <-
        (Rmax[i] * ka * incl_concentrations[i])/(ka*incl_concentrations[i] + kd)

    }

   assoc_formula_second_term <-
      (1 - exp(-((ka*incl_concentrations[i] + kd)*(df_assoc$Time + t0[i]))))

    assoc_formula_full <- assoc_formula_first_term * assoc_formula_second_term

    err_assoc <- c(err_assoc, df_assoc$RU - assoc_formula_full)

    end_of_association_RU <-  assoc_formula_first_term *
      (1 - exp(-((ka*incl_concentrations[i] + kd)*(association + t0[i]))))

    dissoc_decay <- exp(-kd*(df_dissoc$Time - association))

    if (bulkshift)
      end_of_association_RU <- end_of_association_RU + shift[i]

    dissoc_formula_full <- end_of_association_RU * dissoc_decay

    err_dissoc <-   c(err_dissoc, df_dissoc$RU - dissoc_formula_full)

  }
  c(err_assoc, err_dissoc)
}

# Internal. Called by fit_association_dissociation

get_fit_outcomes <- function(Rmax, ka, t0, kd, df, num_conc,
                             incl_concentrations, association, shift, global_rmax){

  full_output_RU <- NULL

  for (i in 1:num_conc){

    df_i <- df %>% dplyr::filter(Concentration == incl_concentrations[i])
 #   RU <- df_i$RU
    Time <- df_i$Time
    Concentration <- df_i$Concentration

    df_i %>% dplyr::filter(AssocIndicator == 1) -> df_assoc
    df_i %>% dplyr::filter(DissocIndicator == 1) -> df_dissoc

    if (global_rmax){
      assoc_formula_first_term <-
        (Rmax * ka * df_assoc$Concentration)/(ka*df_assoc$Concentration + kd)
      end_of_association_RU <-  (Rmax * ka * df_dissoc$Concentration)/(ka*df_dissoc$Concentration + kd) *
        (1 - exp(-((ka*df_dissoc$Concentration + kd)*(association + t0[i]))))


    } else{
      assoc_formula_first_term <-
        (Rmax[i] * ka * df_assoc$Concentration)/(ka*df_assoc$Concentration + kd)
      end_of_association_RU <-  (Rmax[i] * ka * df_dissoc$Concentration)/(ka*df_dissoc$Concentration + kd) *
        (1 - exp(-((ka*df_dissoc$Concentration + kd)*(association + t0[i]))))


    }
     assoc_formula_second_term <-
      (1 - exp(-((ka*df_assoc$Concentration + kd)*(df_assoc$Time + t0[i]))))

    assoc_formula_full <- assoc_formula_first_term * assoc_formula_second_term

    df_assoc$RU <- assoc_formula_full


    dissoc_decay <- exp(-kd*(df_dissoc$Time - association))

    # shift is zero if no bulkshift

    end_of_association_RU <- end_of_association_RU + shift[i]

    dissoc_formula_full <- end_of_association_RU * dissoc_decay

    df_dissoc$RU <- dissoc_formula_full

    full_output_RU <- dplyr::bind_rows(full_output_RU, df_assoc, df_dissoc)

  }
  # return fitted values
  full_output_RU %>% dplyr::select(Time, RU, Concentration)
}

#' Fit sensorgrams for a given well.
#' @param well_idx The corresponding well in the extended sample sheet.
#' @param sample_info A tibble. The expanded sample sheet.
#' @param x_vals A tibble. Time columns from the Carterra output for all concentrations.
#' @param y_vals A tibble RU values from the Carterra output for all concentrations.
#' @param incl_conc_values A vector. The number of concentrations to be included in this plot.
#' @param min_allowed_kd A number. The lowest kd value that can be reliably fit.
#' The default is $10^{-5}$
#' @param max_iterations A number. The maximum number of iterations to perform if the
#' fit does not converge.
#' @param ptol A number. The tolerance level to determine convergence of parameters.
#' @param ftol A number. The tolerance level to determine convergence of the error function.
#' @return A list.
#'
#' fit_association_dissociation <- function(well_idx, sample_info, x_vals, y_vals,
#' incl_concentrations_values, min_allowed_kd = 10^(-5),
#' max_iterations = 500, ptol = 10^(-10), ftol = 10^(-10))


fit_association_dissociation <- function(well_idx, sample_info, x_vals, y_vals,
                            incl_concentrations_values, n_time_points,
                            min_allowed_kd = 10^(-5),
                            max_iterations = 500,
                            ptol = 10^(-10),
                            ftol = 10^(-10)){

  # this function will fit all selected concentrations for one well

  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`

  if (sample_info[well_idx,]$Bulkshift == "Y")
    bulkshift <- TRUE else
    bulkshift <- FALSE

  if (sample_info[well_idx,]$`Global Rmax` == "Y")
    global_rmax <- TRUE else
      global_rmax <- FALSE

  start_idx <- sample_info[well_idx,]$FirstInclConcIdx
  num_conc <- sample_info[well_idx,]$NumInclConc
  end_idx <- start_idx + num_conc - 1

  association <- sample_info[well_idx,]$Association
  dissociation <- sample_info[well_idx,]$Dissociation

  assoc_start <- baseline + baseline_start
  assoc_end <- assoc_start + association

  dissoc_start <- assoc_end
  dissoc_end <- assoc_end + dissociation

  if (sample_info[well_idx,]$`Automate Dissoc. Window` == "Y" & !(is.na(sample_info[well_idx, ]$DissocEnd)))
    dissoc_end <- sample_info[well_idx,]$DissocEnd

  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)

  Time <- x_vals[, start_idx:end_idx] %>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select(value)
  RU <- y_vals[, start_idx:end_idx]%>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select("value")

  incl_concentrations <-
    incl_concentrations_values[start_idx:end_idx]

  purrr::map_dfr(.x = tibble::tibble(incl_concentrations),
          .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    dplyr::arrange(incl_concentrations) -> incl_concentrations_rep

  df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU, "Concentration" = incl_concentrations_rep))
  colnames(df) <- c("Time", "RU", "Concentration") #force correct names - dplyr is doing weird things


  #do both dissociation and association
  df %>% dplyr::mutate(AssocIndicator =
                  ifelse((Time >= assoc_start & Time < assoc_end), 1, 0),
                         DissocIndicator = ifelse(Time > dissoc_start & Time < dissoc_end, 1, 0)) -> df
  df %>% dplyr::filter((AssocIndicator == 1 | DissocIndicator == 1)) -> df

  df %>% dplyr::group_by(Concentration) %>% dplyr::summarise(max = max(RU, na.rm = TRUE)) -> Rmax_start_df
  df %>% dplyr::group_by(Concentration) %>% dplyr::summarise(min = min(RU, na.rm = TRUE)) -> R0_start
  t0_start <- rep(0, num_conc)

  if (global_rmax){
    Rmax_start <- max(Rmax_start_df$max, na.rm = TRUE)
  } else {
    Rmax_start <- Rmax_start_df$max
  }

  kd_start <- 10^(-5)
  ka_start <-  10^(5)

  # shift time to start at zero
  df$Time <- df$Time - assoc_start

  # bulkshift initial value
  shift <- rep(0, num_conc)

  if (bulkshift){
    init_params <- c(Rmax_start, ka_start, t0_start, kd_start, shift)
    fit_result <- minpack.lm::nls.lm(init_params,fn = fit_as_system, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc, association = association,
                         bulkshift,
                         global_rmax = global_rmax,
                         control = minpack.lm::nls.lm.control(maxiter = max_iterations, ptol = ptol, ftol = ftol),
                         lower = c(rep(0, length(Rmax_start)), 10, rep(-Inf,num_conc), min_allowed_kd, rep(-100, num_conc)),
                         jac = NULL,
                       upper = c(rep(400, length(Rmax_start)), 10^7, rep(Inf,num_conc), 1, rep(100, num_conc)))
                       #upper = NULL)
  } else {
    init_params <- c(Rmax_start, ka_start, t0_start, kd_start)

    fit_result <- minpack.lm::nls.lm(init_params,fn = fit_as_system, df = df,
                         incl_concentrations = incl_concentrations, num_conc = num_conc, association = association,
                         bulkshift,
                         global_rmax = global_rmax,
                         control = minpack.lm::nls.lm.control(maxiter = max_iterations, ptol = ptol, ftol = ftol),
                         lower = c(rep(0, length(Rmax_start)), 10, rep(-Inf,num_conc), min_allowed_kd),
                         jac = NULL,
                         upper = c(rep(400, length(Rmax_start)), 10^7, rep(Inf,num_conc), 1))

  }

  pars <- stats::coefficients(fit_result)

  if (global_rmax){
    #pars (global "Rmax","ka", "tstart" one for each concentration)

    Rmax <- pars[1]
    ka <- pars[2]
    t0 <- pars[3:(2 + num_conc)]
    kd <- pars[(3 + num_conc)]

    if (bulkshift)
      shift <- pars[(num_conc + 4):(2*num_conc + 3)]


  } else {
    #pars ("Rmax" one for each concentration,"ka", "tstart" one for each concentration)

    Rmax <- pars[1:num_conc]
    ka <- pars[num_conc+1]
    t0 <- pars[(num_conc + 2):(2*num_conc + 1)]
    kd <- pars[2*num_conc + 2]

    if (bulkshift)
      shift <- pars[(2*num_conc +3):(3*num_conc + 2)]

  }



  fit_outcomes <- get_fit_outcomes(Rmax, ka, t0, kd, df, num_conc,
                                   incl_concentrations, association = association, shift = shift, global_rmax = global_rmax)

   list("FitResult" = fit_result, "FitOutcomes" = fit_outcomes)
}


combine_output <- function(well_idx, fits_list, plot_list_out, rc_list, sample_info){

  if (!is.null(fits_list[[well_idx]]$error))
    return(NULL)

  if (sample_info[well_idx,]$`Global Rmax` == "Y")
    global_rmax <- TRUE else
      global_rmax <- FALSE

  if (sample_info[well_idx,]$`Bulkshift` == "Y")
    bulkshift <- TRUE else
      bulkshift <- FALSE

  num_conc <- sample_info[well_idx,]$NumInclConc

  if (global_rmax)
    Rmax_label <- "Rmax" else
      Rmax_label <- purrr::map_dfr(tibble::tibble(1:num_conc), function(x) paste("Rmax", x))


  R0_label <- purrr::map_dfr(tibble::tibble(1:num_conc), function(x) paste("R_0", x))
  bulkshift_label <- purrr::map_dfr(tibble::tibble(1:num_conc), function(x) paste("Bulkshift", x))


  pars <- stats::coefficients(fits_list[[well_idx]]$result$FitResult)

  if (bulkshift)
    par_names <- purrr::as_vector(purrr::flatten(c(Rmax_label, "ka", R0_label, "kd", bulkshift_label))) else
    par_names <- purrr::as_vector(purrr::flatten(c(Rmax_label, "ka", R0_label, "kd")))

  # result_summary <- summary(fits_list[[well_idx]]$result$FitResult)
  #R's built-in summary method doesn't play nicely when the some of the parameters hit their limiting values (the hessian is singular)
  # I've adapted the function to return NA's for std error when the limits are reached.

  result_summary <- summary_fit_with_constraints(fits_list[[well_idx]]$result$FitResult)
  #summary_fit_with_constraints returns the coefficients table from summary.minpack.lm

  summary_names <- colnames(result_summary)
  result_summary %>% tibble::as_tibble -> par_err_table

  colnames(par_err_table) <- summary_names
  par_err_table <- suppressMessages(dplyr::bind_cols(Names = par_names, par_err_table))

  par_err_table %>% dplyr::filter(!stringr::str_detect(Names,"R_0")) -> par_err_table
  par_err_table %>% dplyr::filter(!stringr::str_detect(Names,"Bulkshift")) -> par_err_table


  par_names <- par_err_table$Names

  par_err_table %>%
    dplyr::filter(Names == "ka" | Names == "kd") %>%
    dplyr::select(Estimate, `Std. Error`)  %>%
    dplyr::mutate(Estimate = format(signif(Estimate, 3),big.mark=",",decimal.mark=".", scientific = TRUE)) %>%
    dplyr::mutate(`Std. Error` = format(signif(`Std. Error`, 3),big.mark=",",decimal.mark=".", scientific = TRUE)) -> kakd_out

  par_err_table %>%
    dplyr::filter(!(Names == "ka" | Names == "kd")) %>%
    dplyr::select(Estimate, `Std. Error`)  %>%
    dplyr::mutate(Estimate = format(round(Estimate,2),big.mark=",",decimal.mark=".", scientific = FALSE))  %>%
    dplyr::mutate(`Std. Error` = format(round(`Std. Error`, 2),big.mark=",",decimal.mark=".", scientific = FALSE)) -> rest_out

  dplyr::bind_rows(rest_out, kakd_out) %>%
    gridExtra::tableGrob(rows = par_names, theme = gridExtra::ttheme_minimal()) -> tb1

  stats::residuals(fits_list[[well_idx]]$result$FitResult) -> RU_resid
  fits_list[[well_idx]]$result$FitOutcomes$Time -> Time_resid
  fits_list[[well_idx]]$result$FitOutcomes$Concentration -> Concentration_resid


  resid_plot <- ggplot2::ggplot(data = tibble::tibble(Residuals = RU_resid, Time = Time_resid, Concentration = forcats::as_factor(Concentration_resid)),
                       ggplot2::aes(x = Time, y = Residuals, color = Concentration)) + ggplot2::geom_point(size = 0.01) +
                          ggplot2::ggtitle(label = "Residuals")

  gridExtra::grid.arrange(plot_list_out[[well_idx]], tb1, resid_plot,
               rc_list[[well_idx]], ncol=2)
}

print_output <- function(well_idx, pages_list, plot_list_out, sample_info){

  if (is.null(pages_list[[well_idx]]$error) & !is.null(pages_list[[well_idx]]$result))
    return(gridExtra::arrangeGrob(pages_list[[well_idx]]$result))

  err_msg <- paste("The following well has an unrecoverable error:",
                     well_idx, "Block ", sample_info$Block, "Row", sample_info$Row)
  err_msg <- paste0(err_msg, sample_info$Column)

  if (is.null(plot_list_out[[well_idx]]))
     return(p1 = gridExtra::arrangeGrob(grid::textGrob(err_msg), plot_list[[well_idx]]))
  else
     return(p1 = gridExtra::arrangeGrob(grid::textGrob(err_msg), plot_list_out[[well_idx]]))

}
get_response_curve <- function(well_idx, sample_info, x_vals, y_vals,
                               all_concentrations_values,
                               incl_concentrations_values,
                               n_time_points){

  start_incl_idx <- sample_info[well_idx,]$FirstInclConcIdx
  start_idx <- sample_info[well_idx,]$FirstConcIdx
  num_conc <- sample_info[well_idx,]$NumConc
  num_incl_conc <- sample_info[well_idx,]$NumInclConc
  end_incl_idx <- start_incl_idx + num_incl_conc - 1
  end_idx <- start_idx + num_conc - 1

  ligand_desc <- sample_info[well_idx,]$Ligand
  baseline <- sample_info[well_idx,]$Baseline
  baseline_start <- sample_info[well_idx,]$`Bsl Start`
  association <- sample_info[well_idx,]$Association

  n_vals <- dim(x_vals)[2]
  names(x_vals) <- as.character(1:n_vals)
  names(y_vals) <- as.character(1:n_vals)


  Time <- x_vals[, start_idx:end_idx] %>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select(value)
  RU <- y_vals[, start_idx:end_idx]%>%
    tidyr::pivot_longer(cols = tidyselect::everything()) %>% dplyr::arrange(as.numeric(name)) %>%
    dplyr::select("value")


  numerical_concentration <- all_concentrations_values[start_idx:end_idx]
  numerical_concentration_incl <-
    incl_concentrations_values[start_incl_idx:end_incl_idx]

  purrr::map_dfr(.x = tibble::tibble(numerical_concentration), .f = function(x, n_time_points) rep(x, n_time_points), n_time_points) %>%
    dplyr::arrange(numerical_concentration) -> numerical_concentration

  Concentrations <- numerical_concentration

  df <- suppressMessages(dplyr::bind_cols("Time" = Time, "RU" = RU, "Concentration" = Concentrations))

  colnames(df) <- c("Time", "RU", "Concentration")


  df %>% dplyr::filter((Time >= baseline + baseline_start + association - 10)
                & (Time <= baseline + baseline_start + association - 5)) %>%
    dplyr::group_by(Concentration) %>%
    dplyr::summarise(AverageRU = mean(RU, na.rm = TRUE)) -> df_RC

  df_RC %>% dplyr::mutate(Included =
                     forcats::as_factor(ifelse(Concentration %in% numerical_concentration_incl,
                                      "Yes", "No"))) -> df_RC


  ggplot2::ggplot(df_RC, ggplot2::aes(x = Concentration,
                    y = AverageRU)) +
    ggplot2::geom_point(ggplot2::aes(color = Included)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle(ligand_desc)
}

get_csv <- function(well_idx, fits_list, sample_info){

  num_conc <- sample_info[well_idx,]$NumInclConc

  if (sample_info[well_idx,]$`Global Rmax` == "Y")
    global_rmax <- TRUE else
      global_rmax <- FALSE

  if (!global_rmax){
    Rmax <- rep(NA, 5)
    Rmax_se <- rep(NA,5)
  } else {
    Rmax <- NA
    Rmax_se <- NA
  }



  Bulkshift <- rep(NA, 5)
  Bulkshift_se <- rep(NA,5)
  R0 <- rep(NA,5)
  R0_se <- rep(NA,5)

  if (is.null(fits_list[well_idx]$FitResult)){
    return(c(Rmax = Rmax, Rmax_se = Rmax_se, ka = NA, ka_se = NA, kd = NA,
             kd_se = NA, Bulkshift = Bulkshift, Bulkshift_se = Bulkshift_se, R0 = R0, R0_se = R0_se))

  }

  pars <- stats::coefficients(fits_list[[well_idx]]$result$FitResult)

  if (sample_info[well_idx, ]$Bulkshift == "Y")
    bulkshift <- TRUE else
    bulkshift <- FALSE



  result_summary <- summary_fit_with_constraints(fits_list[[well_idx]]$result$FitResult)

  # first column of summary is the estimate. Second is standard error

  summary_names <- colnames(result_summary)
  result_summary %>% tibble::as_tibble %>% dplyr::select(Estimate, `Std. Error`) -> par_err_table

  colnames(par_err_table) <- summary_names[1:2]

  if (global_rmax){
    Rmax <- par_err_table[1,]$Estimate
    Rmax_se <- par_err_table[1,]$`Std. Error`
    ka <- par_err_table[2,]$Estimate
    ka_se <- par_err_table[2,]$`Std. Error`
    R0[1:num_conc] <- par_err_table[3:(num_conc + 2), ]$Estimate
    R0_se[1:num_conc] <- par_err_table[3:(num_conc + 2), ]$`Std. Error`
    kd <- par_err_table[2*num_conc + 3,]$Estimate
    kd_se <- par_err_table[2*num_conc + 3,]$`Std. Error`
    if (bulkshift){
      Bulkshift[1:num_conc] <- par_err_table[(2*num_conc + 4):(3*num_conc+3), ]$Estimate
      Bulkshift_se[1:num_conc] <- par_err_table[(2*num_conc + 4):(3*num_conc+3), ]$`Std. Error`
    }


  } else {
    Rmax[1:num_conc] <- par_err_table[1:num_conc,]$Estimate
    Rmax_se[1:num_conc] <- par_err_table[1:num_conc,]$`Std. Error`
    ka <- par_err_table[num_conc + 1,]$Estimate
    ka_se <- par_err_table[num_conc + 1,]$`Std. Error`
    R0[1:num_conc] <- par_err_table[(num_conc+2):(2*num_conc + 1), ]$Estimate
    R0_se[1:num_conc] <- par_err_table[(num_conc+2):(2*num_conc + 1), ]$`Std. Error`
    kd <- par_err_table[2*num_conc + 2,]$Estimate
    kd_se <- par_err_table[2*num_conc + 2,]$`Std. Error`
    if (bulkshift){
      Bulkshift[1:num_conc] <- par_err_table[(2*num_conc + 3):(3*num_conc+2), ]$Estimate
      Bulkshift_se[1:num_conc] <- par_err_table[(2*num_conc + 3):(3*num_conc+2), ]$`Std. Error`
    }

  }

  c(Rmax = Rmax, Rmax_se = Rmax_se, ka = ka, ka_se = ka_se, kd = kd, kd_se = kd_se, Bulkshift = Bulkshift, Bulkshift_se = Bulkshift_se, R0 = R0, R0_se = R0_se)
}
