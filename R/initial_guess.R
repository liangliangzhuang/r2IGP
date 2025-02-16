#' Initial Log-Likelihood Calculation
#'
#' This function computes the log-likelihood for a given set of parameters for different models.
#' The models "M0", "M1", "M2" involve gamma and lambda terms, while models "M3" and "M4" involve
#' only lambda. The function calculates the log-likelihood based on the observed data.
#'
#' @param params A numeric vector containing initial values for the parameters.
#'   The vector should be of length 5 for models "M0", "M1", "M2", or length 3 for models "M3" and "M4".
#' @param model A character string specifying the model type. Must be one of "M0", "M1", "M2", "M3", or "M4".
#' @param type1 A string specifying the type of distribution or model for the lambda calculation.
#' @param data1 A numeric vector containing the observed data for a single product (response variable).
#' @param t A numeric vector of time or other covariates used in the lambda function.
#' @param u A numeric vector of additional covariates for the lambda function (can be NULL for some models).
#'
#' @return A numeric value representing the negative log-likelihood for the given parameters.
#'
#' @export
init.log.likelihood <- function(params = rep(1, 5), model, type1, data1 = sim_dat$diff_Y_t[, 1], t = t, u = NULL) {
  # Model "M0", "M1", "M2" require a parameter for gamma (params[5])
  if (model %in% c("M0", "M1", "M2")) {
    gamma <- params[5] # The gamma parameter

    # Calculate lambda (both Lambda_t and Lambda_u for models M0, M1, M2)
    lambda <- Lambda_fun(type = type1, par0 = params[1:4], t = t, u = u)
    Lambda_t <- lambda$Lambda_t
    Lambda_u <- lambda$Lambda_u

    # Compute log-likelihood for these models
    loglik <- -0.5 * log(2 * pi) - (3 / 2) * log(data1) + log(Lambda_t + Lambda_u) -
      0.5 * (gamma * sqrt(data1) - (Lambda_t + Lambda_u) / sqrt(data1))^2

    # Models "M3" and "M4" only use Lambda_t, with a different parameterization for gamma
  } else if (model %in% c("M3", "M4")) {
    gamma <- params[3] # The gamma parameter for M3 and M4

    # Calculate lambda (only Lambda_t for models M3 and M4)
    lambda <- Lambda_fun(type = type1, par0 = params[1:2], t = t, u = u)
    Lambda_t <- lambda$Lambda_t

    # Compute log-likelihood for these models
    loglik <- -0.5 * log(2 * pi) - (3 / 2) * log(data1) + log(Lambda_t) -
      0.5 * (gamma * sqrt(data1) - Lambda_t / sqrt(data1))^2
  }

  # Return the negative log-likelihood value
  return(-sum(loglik))
}

#' Initial Parameter Guess Using Optimization
#'
#' This function estimates initial parameter values for different models using optimization techniques.
#' The function uses the `optim` function with the L-BFGS-B method, which supports parameter bounds.
#' It calculates the initial parameter estimates for models "M0", "M1", "M2", "M3", and "M4".
#'
#' @param model A character string specifying the model type. Must be one of "M0", "M1", "M2", "M3", or "M4".
#' @param types A vector specifying the type of distribution or model for the lambda calculation.
#' @param data A numeric matrix with columns representing the data for each product.
#' @param t_list A list of vectors, where each vector contains the time or covariates for each product.
#' @param u_list A list of vectors, where each vector contains the additional covariates for each product (can be NULL).
#' @param init_param A numeric vector containing the initial parameter guesses. Default is a vector of 1's of length 3 or 5.
#'
#' @return A numeric vector containing the initial parameter estimates for the given model.
#'
#' @export
init.guess <- function(model, types, data = sim_dat$diff_Y_t, t_list, u_list, init_param = rep(1, 3)) {
  n <- dim(data)[2]
  m <- dim(data)[1]
  len1 <- length(init_param) # Length of the initial parameter vector
  re <- matrix(NA, n, len1) # Initialize matrix to store results for each product

  # Loop through each product to optimize and find initial parameter estimates
  for (i in 1:n) {
    re[i, ] <- tryCatch(
      optim(init_param, init.log.likelihood,
        model = model, type1 = types,
        t = t_list[[i]], u = u_list[[i]], data1 = data[!is.na(data[, i]), i],
        method = "L-BFGS-B", # Use L-BFGS-B method with bounds on parameters
        lower = rep(0, len1), # Set lower bounds for all parameters to 0
        upper = rep(Inf, len1), # Set upper bounds to infinity
        control = list(maxit = 3000)
      )$par, # Max iterations for the optimizer
      error = function(e) {
        return(rep(NA, len1))
      } # If an error occurs, return NA
    )
  }

  # Define initial parameter estimates based on the model type
  if (model == "M0") {
    inti_par <- c(apply(re, 2, mean, na.rm = TRUE)[1:4], mean(re[, 5], na.rm = TRUE), var(re[, 5], na.rm = TRUE))
  } else if (model == "M1") {
    inti_par <- c(apply(re, 2, mean, na.rm = TRUE)[1:4], mean(re[, 5], na.rm = TRUE))
  } else if (model == "M2") {
    inti_par <- c(
      apply(re, 2, mean, na.rm = TRUE)[1:4], mean(re[, 5], na.rm = TRUE),
      var(re[, 5], na.rm = TRUE) / mean(re[, 5], na.rm = TRUE)^2
    )
  } else if (model == "M3") {
    inti_par <- c(apply(re, 2, mean, na.rm = TRUE)[1:2], mean(re[, 3], na.rm = TRUE), var(re[, 3], na.rm = TRUE))
  } else if (model == "M4") {
    inti_par <- c(apply(re, 2, mean, na.rm = TRUE)[1:2], mean(re[, 3], na.rm = TRUE))
  }

  return(inti_par) # Return the calculated initial parameters
}
