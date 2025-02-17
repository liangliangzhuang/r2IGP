#' Reliability Function Calculation
#'
#' This function computes the reliability function (CDF) based on a given model and parameters.
#' It is designed to handle multiple models (M0, M1, M2, M3, M4) and computes the reliability over a specified range of time.
#'
#' @param model A character string indicating the model type. It can be one of "M0", "M1", "M2", "M3", "M4".
#' @param par A numeric vector of parameters used in the model. The length and meaning of the vector depend on the selected model.
#' @param type A character string indicating the type of distribution or function to use (e.g., "Gamma", "Weibull").
#' @param t_ind A numeric value indicating the time at which to evaluate the CDF.
#' @param u_ind A numeric value indicating the additional variable (e.g., uncertainty, user-defined parameter).
#' @param D A numeric value indicating a scale parameter, default is 300.
#'
#' @return A numeric value representing the reliability at time `t_ind` and additional parameters.
#' @importFrom stats pnorm
#' @export
uncon_lifetime_CDF <- function(model, par = em_para_re, type, t_ind, u_ind, D = 300) {

  # Compute Lambda values
  Lambda_re = Lambda_fun(type = type, par0 = par, t = c(0, t_ind), u = c(0, u_ind))
  Lambda_t = Lambda_re$Lambda_t  # Time-dependent Lambda
  Lambda_u = ifelse(length(Lambda_re$Lambda_u) == 0 | is.na(Lambda_re$Lambda_u), 0, Lambda_re$Lambda_u)  # Uncertainty-dependent Lambda

  # Combine Lambda values
  Lambda = Lambda_t + Lambda_u

  if(model %in% c("M0", "M3", "M2")) {
    # Model M2
    if (model == "M2") {
      kappa = as.numeric(par[5])
      sigma2 = as.numeric(par[6])
      Phi_term1 <- (-Lambda / sqrt(D) + kappa * sqrt(D)) / sqrt(1 + D * sigma2 * kappa^2)
      exp_term <- 2 * Lambda * kappa + 2 * Lambda^2 * sigma2 * kappa^2
      term2 = ifelse(exp_term > 500, 10000, exp(exp_term))  # Safe handling for large exponent
      Phi_term2 <- (-Lambda / sqrt(D) - sqrt(D) * kappa - 2 * Lambda * sqrt(D) * sigma2 * kappa^2) / sqrt(1 + D * sigma2 * kappa^2)
    } else {
      # Model M0 or M3
      kappa = as.numeric(par[5])
      sigma2 = as.numeric(par[6])
      if (model == "M3") {
        kappa = as.numeric(par[3])
        sigma2 = as.numeric(par[4])
      }

      # Compute terms for the CDF
      Phi_term1 <- (-Lambda / sqrt(D) + kappa * sqrt(D)) / sqrt(1 + D * sigma2)
      exp_term <- 2 * Lambda * kappa + 2 * Lambda^2 * sigma2
      term2 = ifelse(exp_term > 500, 10000, exp(exp_term))  # Safe handling for large exponent
      Phi_term2 <- (-Lambda / sqrt(D) - sqrt(D) * kappa - 2 * Lambda * sqrt(D) * sigma2) / sqrt(1 + D * sigma2)
    }

    # Compute the CDF using the normal CDF function pnorm
    F_TD <- 1 - pnorm(Phi_term1) - term2 * pnorm(Phi_term2)

  } else if (model %in% c("M1", "M4")) {
    # Model M1 or M4
    if (model == "M1") {
      gamma = as.numeric(par[5])  # M1 dual-scale
    } else {
      gamma = as.numeric(par[3])  # M4 single-scale
    }

    # Compute the CDF using prIG function (Inverse Gamma distribution)
    F_TD <- 1 - prIG(D, Lambda, gamma)
  }

  return(F_TD)
}


#' Compute Reliability over Time
#'
#' This function calculates the reliability function over a range of time points for a given model and parameters.
#' It can generate a plot of the reliability function and compute the Mean Time to Failure (MTTF) based on the calculated CDF.
#'
#' @param par A numeric vector of parameters used in the model.
#' @param model A character string indicating the model type ("M0", "M1", "M2", "M3", "M4").
#' @param t A numeric vector of time points at which to evaluate the reliability.
#' @param u A numeric vector of uncertainty or additional variables, default is NULL.
#' @param type A character string indicating the type of distribution or function to use (e.g., "Gamma", "Weibull").
#' @param D A numeric value indicating a scale parameter, default is 10.
#' @param bin A numeric value indicating the number of bins to use when discretizing the time range, default is 10.
#'
#' @return A list containing:
#' \item{result}{A data frame with the computed reliability for each time point.}
#' \item{MTTF}{The computed Mean Time to Failure (MTTF).}
#' \item{plot}{A ggplot object showing the reliability over time.}
#' @import ggplot2
#' @export
Reliability <- function(par = em_para_re[[1]], model = models, t, u, type = types, D = 10, bin = 10) {

  # Discretize the time range
  t_fine = seq(0, t[length(t)], 1 / bin)

  if (is.null(u)) {
    u_fine = NULL
  } else {
    u_fine = as.numeric(unlist(sapply(1:(length(u) - 1), function(i) {
      seq(u[i], u[i + 1], length.out = bin + 1)[-bin - 1]
    })))
    u_fine <- c(u_fine, u[length(u)])  # Add the last point
  }

  # Compute reliability for each time point
  uncon = numeric(length(t_fine))

  for (i in 1:length(t_fine)) {
    uncon[i] = uncon_lifetime_CDF(model = model, par = par, type = type, D = D, t_ind = t_fine[i], u_ind = u_fine[i])
  }

  RE = data.frame("Time" = t_fine[-1], "Reliability" = 1 - uncon[-1])  # Compute reliability

  # Plot the reliability function
  p1 = ggplot(RE, aes(Time, Reliability)) +
    geom_line() +
    theme_bw() +
    theme(panel.grid = element_blank())

  # Compute MTTF using trapezoidal rule
  mttf <- trapz(t_fine[-1], uncon[-1])

  return(list("result" = RE, "MTTF" = mttf, "plot" = p1))
}
