#' Function to Filter Rows Based on Criteria
#'
#' This function is used to filter out rows that contain infinite values or values greater than a threshold.
#' It checks whether any of the columns are set to FALSE for EM_scen, Bayes_scen, or interval_scen,
#' and also whether the row contains any infinite values or values exceeding 1000.
#'
#' @param row A numeric vector representing a single row of data.
#'
#' @return A logical value: TRUE if the row should be kept, FALSE otherwise.
#'
#' @export
selected_fun <- function(row) {
  # If any of the scenarios (EM_scen, Bayes_scen, interval_scen) are FALSE, retain the row
  if (any(EM_scen == FALSE | Bayes_scen == FALSE | interval_scen == FALSE)) {
    return(FALSE) # Keep the row if any scenario is FALSE
  } else {
    # Remove rows that contain infinite values or values greater than 1000
    return(any(is.infinite(row)) | any(row > 1000))
  }
}


#' Function to Remove Rows Containing Infinite Values or Exceeding Thresholds
#'
#' This function applies the `selected_fun` to each row of the data frame to filter out rows
#' that contain infinite values or exceed a threshold. It is used to clean the data before
#' further analysis.
#'
#' @param df A data frame or matrix where rows are being filtered.
#'
#' @return A data frame or matrix with rows that meet the conditions removed.
#'
#' @export
move_infinite <- function(df) {
  # Apply the selected_fun to each row and retain rows that don't meet the condition
  df[!apply(df, 1, selected_fun), ]
}

#' Function to Save Simulation Results
#'
#' This function processes and saves the results of a simulation, calculating summary statistics
#' such as relative bias (RB), root mean squared error (RMSE), relative RMSE (RRMSE),
#' coverage probability (CP), and confidence interval length (LEN). The results are saved into
#' a list of data frames corresponding to different statistics.
#'
#' @param mttf_re A list containing results from MTTF simulations.
#' @param sim_re A list containing simulation results with keys RB, MSE, re_MSE, etc.
#' @param MTTF Logical, indicating whether to compute MTTF statistics. Default is TRUE.
#' @param col_names A character vector containing column names for the results. Default includes
#'   "alpha_t", "beta_t", "alpha_u", "beta_u", "kappa", and "sigma2".
#'
#' @return A list of data frames containing the calculated statistics.
#'
#' @export
save.result <- function(mttf_re = mttf_re, sim_re = sim_re, MTTF = TRUE, col_names = c("alpha_t", "beta_t", "alpha_u", "beta_u", "kappa", "sigma2")) {
  row_names <- c("EM", "Init", "Bayes")
  row_names2 <- c("Hessian", "Bootstrap", "Bayes")

  # Calculate MTTF statistics if requested
  if (MTTF == TRUE) {
    MTTF <- cbind(
      apply(move_infinite(mttf_re[[1]]), 2, mean, na.rm = TRUE) * 100,
      sqrt(apply(move_infinite(mttf_re[[2]]), 2, mean, na.rm = TRUE)) * 100,
      sqrt(apply(move_infinite(mttf_re[[3]]), 2, mean, na.rm = TRUE)) * 100
    )
    colnames(MTTF) <- c("RB", "RMSE", "RRMSE")
    rownames(MTTF) <- c("EM", "Bayes")
  } else {
    MTTF <- matrix(NA, 2, 3)
  }

  # Create a list of results for each statistic
  sheets <- list(
    "RB" = {
      mat <- matrix(round(apply(move_infinite(sim_re$RB), 2, mean, na.rm = TRUE) * 100, 2), nrow = 3, byrow = TRUE)
      rownames(mat) <- row_names
      colnames(mat) <- col_names
      mat
    },
    "RMSE" = {
      mat <- matrix(round(sqrt(apply(move_infinite(sim_re$MSE), 2, mean, na.rm = TRUE)) * 100, 2), nrow = 3, byrow = TRUE)
      rownames(mat) <- row_names
      colnames(mat) <- col_names
      mat
    },
    "RRMSE" = {
      mat <- matrix(round(sqrt(apply(move_infinite(sim_re$re_MSE), 2, mean, na.rm = TRUE)) * 100, 2), nrow = 3, byrow = TRUE)
      rownames(mat) <- row_names
      colnames(mat) <- col_names
      mat
    },
    "CP" = {
      mat <- matrix(round(apply(move_infinite(sim_re$CP), 2, mean, na.rm = TRUE) * 100, 2), nrow = 3, byrow = TRUE)
      rownames(mat) <- row_names2
      colnames(mat) <- col_names
      mat
    },
    "LEN" = {
      mat <- matrix(round(apply(move_infinite(sim_re$LEN), 2, mean, na.rm = TRUE) * 100, 2), nrow = 3, byrow = TRUE)
      rownames(mat) <- row_names2
      colnames(mat) <- col_names
      mat
    },
    "MTTF" = MTTF
  )

  return(sheets)
}

#' Function to Compute Performance Metrics (RB, MSE, RRMSE, CP, LEN)
#'
#' This function computes several performance metrics comparing the estimated parameters (`final_par`)
#' to the true parameter values (`true_value`). It calculates relative bias (RB), mean squared error (MSE),
#' relative mean squared error (re_MSE), coverage probability (CP), and length of the confidence interval (LEN).
#'
#' @param true_value A numeric vector of true parameter values.
#' @param final_par A numeric vector or matrix containing the estimated parameters.
#' @param interval_re A data frame containing the confidence interval bounds (Lower, Upper) for the parameters.
#'   Default is NULL.
#'
#' @return A list containing the computed values for RB, MSE, re_MSE, CP, and LEN.
#'
#' @export
performence.compare <- function(true_value = para, final_par = c(em_para_re, inti_par, bayes_para_re), interval_re = NULL) {
  # Calculate relative bias (RB)
  RB <- as.numeric((final_par - true_value) / true_value)

  # Calculate mean squared error (MSE)
  MSE <- as.numeric((final_par - true_value)^2)

  # Calculate relative mean squared error (re_MSE)
  re_MSE <- as.numeric((final_par - true_value)^2 / (true_value)^2)

  # If confidence intervals are provided, calculate CP and LEN
  if (!is.null(interval_re)) {
    # Coverage probability (CP)
    CP <- as.numeric(ifelse(interval_re["Lower"] <= true_value & true_value <= interval_re["Upper"], 1, 0))

    # Length of the confidence interval (LEN)
    LEN <- as.numeric(unlist(interval_re["Upper"] - interval_re["Lower"]))
  } else {
    CP <- NA
    LEN <- NA
  }

  # Return the results as a list
  results <- list(RB = RB, MSE = MSE, re_MSE = re_MSE, CP = CP, LEN = LEN)
  return(results)
}
