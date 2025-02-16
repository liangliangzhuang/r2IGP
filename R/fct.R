#' Reparameterized Inverse Gaussian Random Number Generator (rIG)
#'
#' Generates random samples from a reparameterized inverse Gaussian distribution.
#'
#' @param n1 Number of random samples.
#' @param lambda1 Shape parameter of the distribution.
#' @param gamma1 Scale parameter of the distribution.
#' @return A numeric vector of random samples.
#' @importFrom SuppDists rinvGauss
#' @export
rrIG <- function(n1, lambda1, gamma1) {
  re <- SuppDists::rinvGauss(n1, nu = lambda1 / gamma1, lambda = lambda1^2)
  return(re)
}

#' Probability Density Function of the Reparameterized Inverse Gaussian (drIG)
#'
#' Computes the probability density function (PDF) of the reparameterized inverse Gaussian distribution.
#'
#' @param n1 Numeric vector of values to evaluate the PDF at.
#' @param lambda1 Shape parameter.
#' @param gamma1 Scale parameter.
#' @return A numeric vector of PDF values.
#' @importFrom SuppDists dinvGauss
#' @export
drIG <- function(n1, lambda1, gamma1) {
  re <- SuppDists::dinvGauss(n1, nu = lambda1 / gamma1, lambda = lambda1^2)
  return(re)
}


#' Cumulative Distribution Function of the Reparameterized Inverse Gaussian (prIG)
#'
#' Computes the cumulative distribution function (CDF) of the reparameterized inverse Gaussian distribution.
#'
#' @param p Numeric vector of probabilities to evaluate the CDF at.
#' @param lambda1 Shape parameter.
#' @param gamma1 Scale parameter.
#' @return A numeric vector of CDF values.
#' @importFrom SuppDists pinvGauss
#' @export
prIG <- function(p, lambda1, gamma1) {
  re <- SuppDists::pinvGauss(p, nu = lambda1 / gamma1, lambda = lambda1^2)
  return(re)
}

#' Lambda Function Calculation
#'
#' Calculates different Lambda values for the given parameters and type.
#'
#' @param type A character string specifying the type of Lambda function (e.g., "pp").
#' @param par0 A vector of four parameters: alpha_t, beta_t, alpha_u, beta_u.
#' @param t Numeric vector representing time values.
#' @param u Numeric vector representing another variable (e.g., units or cycles).
#' @return A list containing the Lambda_t, Lambda_u, and time values.
#' @export
Lambda_fun <- function(type, par0, t, u) {
  alpha_t <- par0[1]
  beta_t <- par0[2]
  alpha_u <- par0[3]
  beta_u <- par0[4]
  # 使用approx函数进行插值
  # n_all <- length(t) + (length(t) - 1) * bin
  # interpolated_data <- approx(x = 1:length(t), y = t, n = n_all)$y
  Lambda_t <- Lambda_u <- numeric()
  if (type == "pp") { # beta * t ^ alpha
    for (j in 1:(length(t) - 1)) {
      Lambda_t[j] <- beta_t * (t[j + 1]^alpha_t - t[j]^alpha_t)
      Lambda_u[j] <- beta_u * (u[j + 1]^alpha_u - u[j]^alpha_u)
    }
  }
  return(list("Lambda_t" = Lambda_t, "Lambda_u" = Lambda_u, "Time" = cbind(t[-1], u[-1])))
}

#' Cumulative Lambda Calculation
#'
#' Computes the cumulative Lambda values based on input parameters and model type.
#'
#' @param par A vector containing parameters.
#' @param type A string specifying the model type (e.g., "pp").
#' @param t Time vector.
#' @param u Another variable vector (optional).
#' @return A list containing cumulative Lambda values for both time and another variable.
#' @export
Lambda_cum <- function(par, type, t, u) {
  Lambda_re <- Lambda_fun(type = type, par0 = par[1:4], t = t, u = u)
  delta_Lambda0_hat <- Lambda_re$Lambda_t
  delta_Lambda_hat <- Lambda_re$Lambda_u
  Lambda0_hat <- cumsum(delta_Lambda0_hat)
  Lambda_hat <- cumsum(delta_Lambda0_hat)
  return(list(Lambda0_hat = Lambda0_hat, Lambda_hat = Lambda_hat))
}

#' Lambda Function Derivative Calculation
#'
#' Computes the derivatives of Lambda functions with respect to model parameters.
#'
#' @param type The model type (e.g., "M0").
#' @param model The specific model for differentiation.
#' @param par0 A vector of model parameters.
#' @param t Time vector.
#' @param u Another variable (e.g., cycles or units).
#' @return A list containing derivatives of Lambda functions for both time and another variable.
#' @export
Lambda_fun_der <- function(type, model, par0, t, u) {
  alpha_t <- par0[1]
  beta_t <- par0[2]
  alpha_u <- par0[3]
  beta_u <- par0[4]

  m <- length(t) # Define m as the number of time steps

  der_alpha_t <- der_beta_t <- der_alpha_u <- der_beta_u <- numeric()
  for (j in 1:m) {
    # power # beta * t^alpha
    if (model %in% c("M0", "M1", "M2", "M3", "M4")) {
      der_alpha_t[j] <- beta_t * (t[j + 1]^(alpha_t) * log(t[j + 1]) - t[j]^(alpha_t) * log(t[j]))
      der_beta_t[j] <- t[j + 1]^alpha_t - t[j]^alpha_t # beta * (tij^alpha - tij-1^alpha)
    } else {
      der_alpha_t[j] <- der_beta_t[j] <- NA
    }

    if (model %in% c("M0", "M1", "M2")) {
      der_alpha_u[j] <- beta_u * (u[j + 1]^(alpha_u) * log(u[j + 1]) - u[j]^(alpha_u) * log(u[j]))
      der_beta_u[j] <- u[j + 1]^alpha_u - u[j]^alpha_u
    } else {
      der_alpha_u[j] <- der_beta_u[j] <- NA
    }
  }
  return(list(
    "der_alpha_u" = der_alpha_u, "der_beta_u" = der_beta_u,
    "der_alpha_t" = der_alpha_t, "der_beta_t" = der_beta_t
  ))
}


#' Simulate Degradation Data
#'
#' Simulates the degradation path data based on the specified model and parameters.
#'
#' @param model A character string specifying the model type (e.g., "M0", "M1", etc.).
#' @param type The type of Lambda function (e.g., "pp").
#' @param m The number of steps or units.
#' @param n The number of samples or simulations.
#' @param t Time vector.
#' @param u A second variable (optional).
#' @param par A vector of parameters for the model.
#' @param gamma Scale parameter (optional, defaults to 1).
#' @return A list containing the simulated degradation paths for different models.
#' @export
sim.dat.path <- function(model = "M0", type, m, n, t, u, par, gamma = 1) {
  if (model %in% c("M0", "M1", "M2")) {
    Lambda_t <- Lambda_u <- list()
    diff_X_t <- diff_Y_t <- diff_Z_t <- Y_t <- X_t <- Z_t <- matrix(NA, m, n)
    # 设定扩散参数
    if (model == "M2") {
      dispersion <- par[5] * gamma # kappa * gamma_i
    } else if (model == "M1") {
      dispersion <- rep(gamma, n) # kappa * gamma 固定效应(一个数)
    } else if (model == "M0") {
      dispersion <- gamma
    }
    # 计算增量
    for (i in 1:n) {
      Lambda_re <- Lambda_fun(type = type, par0 = par[1:4], t = t[[i]], u = u[[i]])
      Lambda_t[[i]] <- Lambda_re$Lambda_t # t
      Lambda_u[[i]] <- Lambda_re$Lambda_u # u
      for (j in 1:length(Lambda_re$Lambda_u)) {
        diff_X_t[j, i] <- rrIG(1, Lambda_t[[i]][j], dispersion[i]) # X X_i(t)
        diff_Z_t[j, i] <- rrIG(1, Lambda_u[[i]][j], dispersion[i]) # Z Z_i(u)
        diff_Y_t[j, i] <- diff_X_t[j, i] + diff_Z_t[j, i] # Y_i(t,u) = X_i(t) + Z_i(u)
      }
    }
    diff_Y_t <- data.frame(diff_Y_t)
    colnames(diff_Y_t) <- paste("n", 1:n, sep = "")

    Y_t <- apply(diff_Y_t, 2, cumsum)
    colnames(Y_t) <- paste("n", 1:n, sep = "") # 累乘求退化量
    X_t <- apply(diff_X_t, 2, cumsum)
    colnames(X_t) <- paste("n", 1:n, sep = "") # 累乘求退化量
    Z_t <- apply(diff_Z_t, 2, cumsum)
    colnames(Z_t) <- paste("n", 1:n, sep = "") # 累乘求退化量
  } else if (model == "M3") {
    Lambda_t <- list()
    diff_Y_t <- Y_t <- matrix(NA, m, n)
    for (i in 1:n) {
      Lambda_re <- Lambda_fun(type = type, par0 = c(par[1:2], 0, 0), t = t[[i]], u = NULL) # 注意这里的两个0， 只考虑时间尺度下 alpha_t
      Lambda_t[[i]] <- Lambda_re$Lambda_t # t
      for (j in 1:m) {
        diff_Y_t[j, i] <- rrIG(1, Lambda_t[[i]][j], gamma[i]) # X
      }
    }
    Y_t <- apply(diff_Y_t, 2, cumsum) # 累乘求退化量
    X_t <- Z_t <- diff_X_t <- diff_Z_t <- Lambda_u <- NA
  } else if (model == "M4") {
    Lambda_t <- list()
    diff_Y_t <- Y_t <- matrix(NA, m, n)
    for (i in 1:n) {
      Lambda_re <- Lambda_fun(type = type, par0 = c(par[1:2], 0, 0), t = t[[i]], u = NULL) # 注意这里的两个0， 只考虑时间尺度下 alpha_t
      Lambda_t[[i]] <- Lambda_re$Lambda_t # t
      for (j in 1:m) {
        diff_Y_t[j, i] <- rrIG(1, Lambda_t[[i]][j], gamma) # Y gamma(一个数)
      }
    }
    Y_t <- apply(diff_Y_t, 2, cumsum) # 累乘求退化量
    X_t <- Z_t <- diff_X_t <- diff_Z_t <- Lambda_u <- NA
  }

  return(list(
    "Y_t" = Y_t, "X_t" = X_t, "Z_t" = Z_t,
    "diff_Y_t" = diff_Y_t, "diff_X_t" = diff_X_t, "diff_Z_t" = diff_Z_t,
    "Lambda_u" = Lambda_u, "Lambda_t" = Lambda_t
  ))
}
