#' Integrand Function for the Gaussian-Legendre Integral
#'
#' This function defines the integrand for the Gaussian-Legendre integration method.
#' It calculates the product of two parts:
#' 1. A power function of \( u_{ij} \), and
#' 2. An exponential function involving `y_tilde`, `v_tilde`, and `tau_tilde`.
#'
#' @param u_ij Numeric, a sample value for the integration variable.
#' @param y_tilde Numeric, parameter for the exponential function.
#' @param v_tilde Numeric, parameter for the exponential function.
#' @param tau_tilde Numeric, parameter for the exponential function.
#'
#' @return The value of the integrand at the given `u_ij`.
#'
#' @export
Integrand_fun <- function(u_ij, y_tilde, v_tilde, tau_tilde) {
  # Define the integrand: a combination of a power law and an exponential function
  (1 - u_ij)^(-3 / 2) * u_ij^(-3 / 2) * exp(
    -0.5 / y_tilde * ((v_tilde^2 / (1 - u_ij)) + (tau_tilde^2 / u_ij))
  )
}

#' Monte Carlo Integration
#'
#' This function estimates an integral using the Monte Carlo method. It randomly samples
#' `u_ij` values from a uniform distribution and computes the mean of the integrand
#' values for these samples.
#'
#' @param period Numeric vector, the integration interval, default is `c(0, 1)`.
#' @param n_samples Integer, number of samples to use for the Monte Carlo integration.
#' @param y_tilde Numeric, parameter for the exponential function in the integrand.
#' @param v_tilde Numeric, parameter for the exponential function in the integrand.
#' @param tau_tilde Numeric, parameter for the exponential function in the integrand.
#'
#' @return A list containing:
#'   - `integral_estimate`: Estimated integral using Monte Carlo.
#'   - `time_avg`: Time taken to perform the calculation.
#'
#' @export
mc_integral <- function(period = c(0, 1), n_samples, y_tilde, v_tilde, tau_tilde) {
  t1 <- Sys.time()
  # Step 1: 从 [0,1] 均匀采样 n_samples 个 u_ij
  u_samples <- runif(n_samples, min = period[1], max = period[2])
  # Step 2: 计算被积函数的值
  integrand_values <- sapply(u_samples, Integrand_fun,
    y_tilde = y_tilde,
    v_tilde = v_tilde, tau_tilde = tau_tilde
  )
  # Step 3: 求样本平均
  if (length(y_tilde) == 1) {
    integral_estimate <- mean(integrand_values, na.rm = TRUE)
  } else {
    integral_estimate <- apply(integrand_values, 1, mean, na.rm = TRUE)
  }

  time_avg <- as.numeric(Sys.time() - t1, units = "secs")

  return(list(integral_estimate, time_avg))
}


#' Trapezoidal Rule Integration
#'
#' This function estimates the integral using the trapezoidal rule. It creates a set of evenly spaced
#' `u_ij` values and computes the average of the integrand at the boundaries and at the intermediate points.
#'
#' @param period Numeric vector, the integration interval, default is `c(0, 1)`.
#' @param y_tilde Numeric, parameter for the exponential function in the integrand.
#' @param v_tilde Numeric, parameter for the exponential function in the integrand.
#' @param tau_tilde Numeric, parameter for the exponential function in the integrand.
#' @param n_intervals Integer, the number of intervals for the trapezoidal rule.
#' @param epsilon Numeric, small value to avoid sampling at the boundaries (0 or 1).
#'
#' @return A list containing:
#'   - `integral_estimate`: Estimated integral using the trapezoidal rule.
#'   - `time_avg`: Time taken to perform the calculation.
#'
#' @export
trapezoidal_integral <- function(period = c(0, 1), y_tilde, v_tilde, tau_tilde, n_intervals = 1000, epsilon = 1e-6) {
  t1 <- Sys.time()
  # 创建 u_ij 的均匀分布节点，但避免 u_ij 等于 0 或 1
  u_vals <- seq(period[1] + epsilon, period[2] - epsilon, length.out = n_intervals)
  # 计算被积函数
  integrand_values <- sapply(u_vals, Integrand_fun,
    y_tilde = y_tilde,
    v_tilde = v_tilde, tau_tilde = tau_tilde
  )
  # 使用梯形法计算积分
  if (length(y_tilde) == 1) {
    integral_estimate <- sum((integrand_values[-1] + integrand_values[-length(integrand_values)]) / 2 * diff(u_vals), na.rm = TRUE)
  } else {
    integral_estimate <- apply((integrand_values[, -1] + integrand_values[, -dim(integrand_values)[2]]) / 2 * diff(u_vals), 1, sum, na.rm = TRUE)
  }
  time_avg <- as.numeric(Sys.time() - t1, units = "secs")
  return(list(integral_estimate, time_avg))
}

#' Gaussian-Legendre Integration
#'
#' This function estimates the integral using the Gaussian-Legendre quadrature method.
#' It uses Legendre polynomial roots and weights to compute the integral over the given interval.
#'
#' @param period Numeric vector, the integration interval, default is `c(0, 1)`.
#' @param y_tilde Numeric, parameter for the exponential function in the integrand.
#' @param v_tilde Numeric, parameter for the exponential function in the integrand.
#' @param tau_tilde Numeric, parameter for the exponential function in the integrand.
#' @param n_points Integer, the number of Gaussian points to use for the integration.
#'
#' @return A list containing:
#'   - `integral_estimate`: Estimated integral using the Gaussian-Legendre quadrature.
#'   - `time_avg`: Time taken to perform the calculation.
#'
#' @export
gaussian_legendre_integral <- function(period = c(0, 1), y_tilde, v_tilde, tau_tilde, n_points = 10) {
  t1 <- Sys.time()
  a <- period[1]
  b <- period[2]
  # 获取勒让德多项式的节点和权重
  gauss_points <- gauss.quad(n_points, kind = "legendre")
  # 转换积分区间为 [a, b]，因为高斯积分默认在 [-1, 1]
  u_vals <- (b - a) / 2 * gauss_points$nodes + (b + a) / 2
  weights <- (b - a) / 2 * gauss_points$weights
  # 计算被积函数的值
  integrand_values <- sapply(u_vals, Integrand_fun,
    y_tilde = y_tilde,
    v_tilde = v_tilde, tau_tilde = tau_tilde
  )
  # 使用高斯积分法计算积分
  if (length(y_tilde) == 1) {
    integral_estimate <- sum(integrand_values * weights)
  } else {
    integral_estimate <- integrand_values %*% weights
  }

  time_avg <- as.numeric(Sys.time() - t1, units = "secs")
  return(list(integral_estimate, time_avg))
}

#' Approximate Integration Method Selector
#'
#' This function selects the appropriate integration method (Monte Carlo, Trapezoidal, or Gaussian-Legendre)
#' based on the user's input and calls the corresponding function to compute the integral.
#'
#' @param method Character, the integration method to use ("mc", "ti", or "gl").
#' @param period Numeric vector, the integration interval.
#' @param y_tilde Numeric, parameter for the exponential function in the integrand.
#' @param v_tilde Numeric, parameter for the exponential function in the integrand.
#' @param tau_tilde Numeric, parameter for the exponential function in the integrand.
#' @param n_samples Integer, number of samples for Monte Carlo integration (if `method` is "mc").
#' @param n_intervals Integer, number of intervals for Trapezoidal integration (if `method` is "ti").
#' @param n_points Integer, number of points for Gaussian-Legendre integration (if `method` is "gl").
#' @param epsilon Numeric, small value for Trapezoidal integration (if `method` is "ti").
#'
#' @return The result from the selected integration method.
#'
#' @export
Approx.integration <- function(method, period, y_tilde, v_tilde, tau_tilde, n_samples = NULL, n_intervals = NULL, n_points = NULL, epsilon = NULL) {
  if (method == "mc") {
    # 调用蒙特卡洛积分
    result <- mc_integral(period = period, n_samples = n_samples, y_tilde = y_tilde, v_tilde = v_tilde, tau_tilde = tau_tilde)
    return(result)
  } else if (method == "ti") {
    # 调用梯形积分
    result <- trapezoidal_integral(period = period, y_tilde = y_tilde, v_tilde = v_tilde, tau_tilde = tau_tilde, n_intervals = n_intervals, epsilon = epsilon)
    return(result)
  } else if (method == "gl") {
    # 调用高斯-勒让德积分
    result <- gaussian_legendre_integral(period = period, y_tilde = y_tilde, v_tilde = v_tilde, tau_tilde = tau_tilde, n_points = n_points)
    return(result)
  }
}


#' Estimate Function for f(u_ij | D)
#'
#' This function computes the ratio \( f(u_{ij} | D) \), where \( f(u_{ij} | D) \) is the
#' likelihood of a given value of \( u_{ij} \) conditioned on the observed data.
#'
#' @param u_ij Numeric, the integration variable.
#' @param method Character, the integration method to use ("mc", "ti", "gl").
#' @param period Numeric vector, the integration period.
#' @param y_tilde Numeric, parameter for the exponential function in the integrand.
#' @param v_tilde Numeric, parameter for the exponential function in the integrand.
#' @param tau_tilde Numeric, parameter for the exponential function in the integrand.
#' @param n_samples Integer, number of samples for Monte Carlo integration.
#' @param n_intervals Integer, number of intervals for Trapezoidal integration.
#' @param n_points Integer, number of points for Gaussian-Legendre integration.
#'
#' @return The likelihood estimate \( f(u_{ij} | D) \).
#'
#' @export
f_uz_given_D_gaussian <- function(u_ij, method = "gl", period = c(0, 1), y_tilde, v_tilde, tau_tilde,
                                  n_samples = NULL, n_intervals = NULL, n_points = 30) {
  # 根据不同的method，可以设置参数，默认为gl。例如：
  # f_uz_given_D_gaussian(0.3, method = "gl", y_tilde = y_tilde, v_tilde = v_tilde,
  #                                         tau_tilde = tau_tilde, n_points = n_points)
  # f_uz_given_D_gaussian(0.3, method = "mc", y_tilde = y_tilde, v_tilde = v_tilde,
  #                       tau_tilde = tau_tilde, n_samples = mc_sample)
  #
  # f_uz_given_D_gaussian(0.3, method = "ti", y_tilde = y_tilde, v_tilde = v_tilde,
  #                       tau_tilde = tau_tilde, n_intervals = ti_intervals)
  meth <- times <- as.numeric()
  # 计算分子
  numerator <- Integrand_fun(u_ij, y_tilde = y_tilde, v_tilde = v_tilde, tau_tilde = tau_tilde)
  # 计算分母
  app_re <- Approx.integration(
    method = method, period = period, y_tilde = y_tilde, v_tilde = v_tilde, tau_tilde = tau_tilde,
    n_samples = n_samples, n_intervals = n_intervals, n_points = n_points, epsilon = 1e-6
  )
  # 返回最终结果
  result <- numerator / app_re[[1]]

  return(result)
}

#' EM Algorithm: Expectation Step (E-step)
#'
#' This function calculates the posterior distribution of the latent variables and computes the expected values
#' for different model configurations. It approximates expectations using methods like Gaussian-Legendre (GL),
#' Monte Carlo (MC), or Trapezoidal (TI) methods.
#'
#' @param serial Integer, specifies the expected value to calculate (e.g., 1 to 6).
#' @param model Character, specifies the model type ("M0", "M1", "M2", etc.).
#' @param method Character, specifies the numerical method for approximation ("gl", "mc", or "ti").
#' @param period Numeric vector of length 2, defines the integration limits (e.g., `c(0, 1)`).
#' @param par1 List, model parameters for each variable. Includes `delta_Lambda` and `delta_Lambda0`.
#' @param gamma_par Numeric vector, parameters for the gamma distribution, used for specific models (e.g., `c(2, 0.1)`).
#' @param data Matrix, observed data where rows represent samples and columns represent variables.
#' @param n_samples Integer, the number of random samples for Monte Carlo integration. Default is NULL.
#' @param n_intervals Integer, the number of intervals for Trapezoidal integration. Default is NULL.
#' @param n_points Integer, the number of points for Gaussian-Legendre integration. Default is NULL.
#'
#' @return Numeric matrix, expected values for the latent variables based on the provided `serial` and `model`.
#'         The result varies depending on the type of expectation (e.g., posterior mean, variance, etc.).
#'
#' @export
E_z <- function(serial, model, method = "gl", period = c(0, 1), par1, gamma_par = c(2, 0.1), data,
                n_samples = NULL, n_intervals = NULL, n_points = NULL) {
  # 计算隐变量的后验分布 p(Z|Y); par = c(alpha0, beta0, alpha, beta, gamma)
  # 求条件密度函数
  n <- dim(data)[2]
  m <- dim(data)[1]
  # 不同方法近似 =======
  ## 1. 设置节点和权重
  if (method == "gl") {
    a <- period[1]
    b <- period[2]
    gauss_points <- gauss.quad(n_points, kind = "legendre")
    # 转换积分区间为 [a, b]，因为高斯积分默认在 [-1, 1]
    u_vals <- (b - a) / 2 * gauss_points$nodes + (b + a) / 2
    weights <- (b - a) / 2 * gauss_points$weights
  } else if (method == "mc") {
    u_vals <- runif(n_samples, min = period[1], max = period[2])
  } else if (method == "ti") {
    u_vals <- seq(period[1] + 1e-6, period[2] - 1e-6, length.out = n_intervals)
  }
  ## 2. 计算不同情景的期望（4 + 2种）
  E_re <- matrix(NA, m, n)
  if (serial %in% 1:4) {
    for (i in 1:n) {
      integrand_values <- sapply(u_vals, function(u_ij) {
        # 根据 serial 参数选择计算逻辑
        if (serial == 1) {
          return(log(1 - u_ij) * f_uz_given_D_gaussian(u_ij,
            method = method, y_tilde = data[, i],
            v_tilde = par1[[i]]$delta_Lambda, tau_tilde = par1[[i]]$delta_Lambda0,
            n_points = n_points, n_samples = n_samples, n_intervals = n_intervals
          ))
        } else if (serial == 2) {
          return((1 / (1 - u_ij)) * f_uz_given_D_gaussian(u_ij,
            method = method, y_tilde = data[, i],
            v_tilde = par1[[i]]$delta_Lambda, tau_tilde = par1[[i]]$delta_Lambda0,
            n_points = n_points, n_samples = n_samples, n_intervals = n_intervals
          ))
        } else if (serial == 3) {
          return((1 / u_ij) * f_uz_given_D_gaussian(u_ij,
            method = method, y_tilde = data[, i],
            v_tilde = par1[[i]]$delta_Lambda, tau_tilde = par1[[i]]$delta_Lambda0,
            n_points = n_points, n_samples = n_samples, n_intervals = n_intervals
          ))
        } else if (serial == 4) {
          return(log(u_ij) * f_uz_given_D_gaussian(u_ij,
            method = method, y_tilde = data[, i],
            v_tilde = par1[[i]]$delta_Lambda, tau_tilde = par1[[i]]$delta_Lambda0,
            n_points = n_points, n_samples = n_samples, n_intervals = n_intervals
          ))
        }
      })

      # 根据 method 进行积分计算
      if (method == "gl") {
        E_re[, i] <- integrand_values %*% weights # 高斯-勒让德积分
      } else if (method == "mc") {
        E_re[, i] <- apply(integrand_values, 1, mean, na.rm = TRUE) # 蒙特卡洛积分
      } else if (method == "ti") {
        E_re[, i] <- apply((integrand_values[, -1] + integrand_values[, -dim(integrand_values)[2]]) / 2 * diff(u_vals), 1, sum, na.rm = TRUE) # 梯形法
      }
    }
  } else if (serial == 5) {
    # 期望 5： E[gamma_i]
    E_re <- as.numeric()
    if (model == "M0") {
      for (i in 1:n) { # 期望
        numerator <- gamma_par[1] / gamma_par[2] + sum(par1[[i]]$delta_Lambda + par1[[i]]$delta_Lambda0)
        denominator <- 1 / gamma_par[2] + sum(data[, i], na.rm = T)
        E_re[i] <- numerator / denominator
      }
    } else if (model == "M2") {
      for (i in 1:n) { # 期望
        numerator <- gamma_par[1] * sum(par1[[i]]$delta_Lambda + par1[[i]]$delta_Lambda0) + 1 / gamma_par[2]
        denominator <- gamma_par[1]^2 * sum(data[, i], na.rm = T) + 1 / gamma_par[2]
        E_re[i] <- numerator / denominator
      }
    } else if (model == "M3") {
      for (i in 1:n) { # 期望
        numerator <- gamma_par[1] / gamma_par[2] + sum(par1[[i]]$delta_Lambda0)
        denominator <- 1 / gamma_par[2] + sum(data[, i], na.rm = T)
        E_re[i] <- numerator / denominator
      }
    } else {
      E_re <- NA
    }
  } else if (serial == 6) {
    # 期望 6： E[gamma^2_i]
    E_re <- as.numeric()
    if (model == "M0") {
      for (i in 1:n) {
        numerator <- 1 / gamma_par[2] + sum(data[, i], na.rm = T) + (gamma_par[1] / gamma_par[2] + sum(par1[[i]]$delta_Lambda + par1[[i]]$delta_Lambda0))^2
        denominator <- (1 / gamma_par[2] + sum(data[, i], na.rm = T))^2
        E_re[i] <- numerator / denominator
      }
    } else if (model == "M2") {
      for (i in 1:n) { # 期望
        Mean_re <- (gamma_par[1] * sum(par1[[i]]$delta_Lambda + par1[[i]]$delta_Lambda0) + 1 / gamma_par[2]) / (gamma_par[1]^2 * sum(data[, i], na.rm = T) + 1 / gamma_par[2])
        Var_re <- 1 / (gamma_par[1]^2 * sum(data[, i], na.rm = T) + 1 / gamma_par[2])
        E_re[i] <- Mean_re^2 + Var_re
      }
    } else if (model == "M3") {
      for (i in 1:n) { # 期望
        numerator <- 1 / gamma_par[2] + sum(data[, i], na.rm = T) + (gamma_par[1] / gamma_par[2] + sum(par1[[i]]$delta_Lambda0))^2
        denominator <- (1 / gamma_par[2] + sum(data[, i], na.rm = T))^2
        E_re[i] <- numerator / denominator
      }
    } else {
      E_re <- NA
    }
  }
  return(E_re)
}

#' EM Algorithm: Expectation-Maximization (EM) Algorithm
#'
#' This function implements the Expectation-Maximization (EM) algorithm, which iteratively optimizes model parameters
#' by alternating between an E-step (Expectation step) and an M-step (Maximization step). The algorithm is used
#' for parameter estimation in latent variable models, with different model types like "M0", "M1", "M2", and "M3".
#'
#' The E-step calculates the expected values of latent variables, and the M-step updates the model parameters
#' based on these expectations. The function also supports different numerical methods for approximating the
#' expectations (e.g., Monte Carlo, Trapezoidal, and Gaussian-Legendre methods).
#'
#' @param data Numeric matrix of size `m x n`, where `m` is the number of samples and `n` is the number of variables.
#' @param par0 Numeric vector of initial model parameters.
#' @param types Character vector specifying the model types for each variable.
#' @param tol1 Numeric value, the tolerance for convergence. The algorithm stops if the parameter change is below this threshold.
#' @param max_iter Integer, the maximum number of iterations for the EM algorithm.
#' @param t_list List, contains the `t` values for each variable (used in Lambda calculations).
#' @param u_list List, contains the `u` values for each variable (used in Lambda calculations).
#' @param model Character, the model type ("M0", "M1", "M2", "M3", etc.) to be used in the estimation.
#' @param n_points Integer, the number of points for Gaussian-Legendre integration.
#' @param n_samples Integer, the number of random samples for Monte Carlo integration.
#' @param n_intervals Integer, the number of intervals for Trapezoidal integration.
#' @param approx.method Character, specifies the numerical method for approximating the expectation ("gl" for Gaussian-Legendre, "mc" for Monte Carlo, or "ti" for Trapezoidal).
#' @param period Numeric vector of length 2, defining the integration limits (e.g., `c(0,1)`).
#' @param show_progress Logical, if TRUE, shows a progress bar during iterations.
#'
#' @return A list containing:
#'   - `par_re`: A matrix of estimated parameters after the final iteration.
#'   - `iter`: The number of iterations performed.
#'   - `mse_all`: The mean squared error between parameter estimates in consecutive iterations, used as a convergence criterion.
#'
#' @export
#' @import utils
#' @importFrom  stats uniroot
#' @importFrom statmod gauss.quad
EM <- function(data, par0, types, tol1 = 0.001, max_iter, t_list, u_list, model,
               n_points = n_points, n_samples = 100, n_intervals = NULL, approx.method,
               period = c(0, 1), show_progress = TRUE) {
  # 根据 `show_progress` 参数来决定是否显示进度条
  if (show_progress) {
    pb <- txtProgressBar(min = 1, max = max_iter, style = 3)
  }
  n <- dim(data)[2]
  m <- dim(data)[1]
  # 设置初始参数
  iter <- 1
  par_re <- matrix(NA, max_iter, length(par0))
  par_re[1, ] <- par0
  par <- list()
  while (iter < max_iter) {
    # 迭代次数
    iter <- iter + 1
    # 参数设置
    init_para <- par_re[iter - 1, ]
    for (i in 1:n) {
      Lambda_re <- Lambda_fun(type = types, par0 = init_para, t = t_list[[i]], u = u_list[[i]])
      par[[i]] <- list("delta_Lambda0" = Lambda_re$Lambda_t, "delta_Lambda" = Lambda_re$Lambda_u)
    }

    # E-step: 计算隐变量的后验分布 p(Z|Y) + 对应期望============
    if (model %in% c("M0", "M1", "M2")) {
      E2 <- E_z(
        serial = 2, model = model, method = approx.method, period = period, par1 = par, data = data,
        n_samples = n_samples, n_intervals = n_intervals, n_points = n_points
      )
      E3 <- E_z(
        serial = 3, model = model, method = approx.method, period = period, par1 = par, data = data,
        n_samples = n_samples, n_intervals = n_intervals, n_points = n_points
      )
    } else {
      E2 <- E3 <- matrix(NA, m, n)
    }

    if (model %in% c("M0", "M2", "M3")) {
      if (model %in% c("M0", "M2")) {
        E5 <- E_z(
          serial = 5, model = model, method = approx.method, period = period, par1 = par, data = data,
          n_samples = n_samples, n_intervals = n_intervals, n_points = n_points,
          gamma_par = init_para[5:6]
        )
        E6 <- E_z(
          serial = 6, model = model, method = approx.method, period = period, par1 = par, data = data,
          n_samples = n_samples, n_intervals = n_intervals, n_points = n_points,
          gamma_par = init_para[5:6]
        )
      } else if (model == "M3") {
        E5 <- E_z(
          serial = 5, model = model, method = approx.method, period = period, par1 = par, data = data,
          n_samples = n_samples, n_intervals = n_intervals, n_points = n_points,
          gamma_par = init_para[3:4]
        )
        E6 <- E_z(
          serial = 6, model = model, method = approx.method, period = period, par1 = par, data = data,
          n_samples = n_samples, n_intervals = n_intervals, n_points = n_points,
          gamma_par = init_para[3:4]
        )
      }
    } else {
      E5 <- E6 <- rep(NA, n)
    }

    # M-step: ===============
    init_para_new <- init_para # 估计参数替换上去
    # 1. alpha_t ====
    alpha_t_prime <- function(x, raw_par = init_para_new, model, data, E = list(E3, E5), types, t_list, u_list) {
      # para_est: beta * t^alpha； 这个里面的power上的alpha
      # 重新修改参数
      raw_par[1] <- x
      part2_re <- ree <- numeric()
      nu_hat <- tau_hat <- lam_der <- list()
      # 重新计算结果（Lambda）
      for (i in 1:n) {
        Lambda_re <- Lambda_fun(type = types, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])
        nu_hat[[i]] <- Lambda_re$Lambda_u
        tau_hat[[i]] <- Lambda_re$Lambda_t
        lam_der[[i]] <- Lambda_fun_der(type = types, model = model, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])$der_alpha_t # beta * t^alpha

        if (model == "M0") {
          part2_re <- lam_der[[i]] * (E[[2]][i] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]])
        } else if (model == "M1") {
          part2_re <- lam_der[[i]] * (raw_par[5] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]]) # raw_par[5] -> gamma
        } else if (model == "M2") {
          part2_re <- lam_der[[i]] * (raw_par[5] * E[[2]][i] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]]) # raw_par[5] -> kappa
        } else if (model == "M3") {
          part2_re <- lam_der[[i]] * (E[[2]][i] - tau_hat[[i]] / data[, i] + 1 / tau_hat[[i]]) # raw_par[5] -> kappa
        } else {
          part2_re <- 0
        }

        ree[i] <- sum(part2_re, na.rm = T)
      }
      return(sum(ree, na.rm = T))
    }
    # 2. beta_t ====
    beta_t_prime <- function(x, raw_par = init_para_new, model, data, E = list(E3, E5), types, t_list, u_list) {
      # para_est: beta * t^alpha； 这个里面的power上的beta
      # 重新修改参数
      raw_par[2] <- x
      part2_re <- ree <- numeric()
      nu_hat <- tau_hat <- lam_der <- list()
      # 重新计算结果（Lambda）
      for (i in 1:n) {
        Lambda_re <- Lambda_fun(type = types, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])
        nu_hat[[i]] <- Lambda_re$Lambda_u
        tau_hat[[i]] <- Lambda_re$Lambda_t
        lam_der[[i]] <- Lambda_fun_der(type = types, model = model, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])$der_beta_t # beta * t^alpha

        if (model == "M0") {
          part2_re <- lam_der[[i]] * (E[[2]][i] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]])
        } else if (model == "M1") {
          part2_re <- lam_der[[i]] * (raw_par[5] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]]) # raw_par[5] -> gamma
        } else if (model == "M2") {
          part2_re <- lam_der[[i]] * (raw_par[5] * E[[2]][i] - tau_hat[[i]] / data[, i] * E[[1]][, i] + 1 / tau_hat[[i]]) # raw_par[5] -> kappa
        } else if (model == "M3") {
          part2_re <- lam_der[[i]] * (E[[2]][i] - tau_hat[[i]] / data[, i] + 1 / tau_hat[[i]]) # raw_par[5] -> kappa
        } else {
          part2_re <- 0
        }

        ree[i] <- sum(part2_re, na.rm = T)
      }
      return(sum(ree, na.rm = T))
    }
    # 3. alpha_u ===
    alpha_u_prime <- function(x, raw_par = init_para_new, data, model, E = list(E2, E5), types, t_list, u_list) {
      # para_est: beta * t^alpha； 这个里面的power上的alpha
      # 重新修改参数
      if (model %in% c("M0", "M1", "M2")) {
        raw_par[3] <- x
        part2_re <- ree <- numeric()
        nu_hat <- tau_hat <- lam_der <- list()
        # 重新计算结果（Lambda）
        for (i in 1:n) {
          Lambda_re <- Lambda_fun(type = types, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])
          nu_hat[[i]] <- Lambda_re$Lambda_u
          tau_hat[[i]] <- Lambda_re$Lambda_t
          lam_der[[i]] <- Lambda_fun_der(type = types, model = model, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])$der_alpha_u # beta * t^alpha

          if (model == "M0") {
            part2_re <- lam_der[[i]] * (E[[2]][i] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]])
          } else if (model == "M1") {
            part2_re <- lam_der[[i]] * (raw_par[5] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]]) # raw_par[5] -> gamma
          } else if (model == "M2") {
            part2_re <- lam_der[[i]] * (raw_par[5] * E[[2]][i] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]]) # raw_par[5] -> kappa
          }
          ree[i] <- sum(part2_re, na.rm = T)
        }
        return(sum(ree, na.rm = T))
      } else {
        return(NA)
      }
    }
    # 4. beta_u ====
    beta_u_prime <- function(x, raw_par = init_para_new, data, model, E = list(E2, E5), types, t_list, u_list) {
      if (model %in% c("M0", "M1", "M2")) {
        # 重新修改参数
        raw_par[4] <- x
        part2_re <- ree <- numeric()
        nu_hat <- tau_hat <- lam_der <- list()
        # 重新计算结果（Lambda）
        for (i in 1:n) {
          Lambda_re <- Lambda_fun(type = types, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])
          nu_hat[[i]] <- Lambda_re$Lambda_u
          tau_hat[[i]] <- Lambda_re$Lambda_t
          lam_der[[i]] <- Lambda_fun_der(type = types, model = model, par0 = raw_par, t = t_list[[i]], u = u_list[[i]])$der_beta_u # beta * t^alpha
          if (model == "M0") {
            part2_re <- lam_der[[i]] * (E[[2]][i] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]])
          } else if (model == "M1") {
            part2_re <- lam_der[[i]] * (raw_par[5] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]]) # raw_par[5] -> gamma
          } else if (model == "M2") {
            part2_re <- lam_der[[i]] * (raw_par[5] * E[[2]][i] - nu_hat[[i]] / data[, i] * E[[1]][, i] + 1 / nu_hat[[i]]) # raw_par[5] -> kappa
          }
          ree[i] <- sum(part2_re, na.rm = T)
        }
        return(sum(ree, na.rm = T))
      } else {
        return(NA)
      }
    }

    # 1. alpha_t计算求导 ====
    hat_alpha_t <- tryCatch(
      uniroot(alpha_t_prime,
        raw_par = init_para_new, data = data, model = model,
        types = types, t_list = t_list, u_list = u_list, E = list(E3, E5), c(0.001, 4)
      )$root,
      error = function(e) {
        return(NA)
      }
    )
    init_para_new[1] <- hat_alpha_t # 修改结果

    # 2. beta_t 计算求导 ====
    hat_beta_t <- tryCatch(uniroot(beta_t_prime, raw_par = init_para_new, data = data, model = model, types = types, t_list = t_list, u_list = u_list, c(0.001, 10))$root,
      error = function(e) {
        return(NA)
      }
    )
    init_para_new[2] <- hat_beta_t # 修改结果


    if (model %in% c("M0", "M1", "M2")) {
      # 3. alpha_u 计算求导 ====
      hat_alpha_u <- tryCatch(uniroot(alpha_u_prime, raw_par = init_para_new, data = data, model = model, types = types, t_list = t_list, u_list = u_list, c(0.001, 4))$root,
        error = function(e) {
          return(NA)
        }
      )
      init_para_new[3] <- hat_alpha_u # 修改结果

      # 4. beta_u 计算求导 ====
      hat_beta_u <- tryCatch(uniroot(beta_u_prime, raw_par = init_para_new, data = data, model = model, types = types, t_list = t_list, u_list = u_list, c(0.001, 10))$root,
        error = function(e) {
          return(NA)
        }
      )

      init_para_new[4] <- hat_beta_u # 修改结果
    }

    # 5. kappa, sigma2 ====
    if (model %in% c("M0", "M2", "M3")) {
      if (model == "M0") {
        hat_kappa <- mean(E5)
        hat_sigma2 <- ifelse(mean(E6 - 2 * hat_kappa * E5) + hat_kappa^2 > 0, mean(E6 - 2 * hat_kappa * E5) + hat_kappa^2, 0)
        init_para_new[5:6] <- c(hat_kappa, hat_sigma2) # 修改结果
      } else if (model == "M2") {
        nu_hat1 <- tau_hat1 <- list()
        term1 <- as.numeric()
        for (i in 1:n) {
          Lambda_re <- Lambda_fun(type = types, par0 = init_para_new, t = t_list[[i]], u = u_list[[i]])
          nu_hat1[[i]] <- Lambda_re$Lambda_u
          tau_hat1[[i]] <- Lambda_re$Lambda_t
          term1[i] <- sum(nu_hat1[[i]] + tau_hat1[[i]])
        }
        term1 <- sum(term1 * E5)
        term2 <- sum(apply(data, 2, sum, na.rm = T) * E6)
        hat_kappa <- term1 / term2
        # hat_kappa = 5
        hat_sigma2 <- mean(E6 - 2 * E5) + 1
        init_para_new[5:6] <- c(hat_kappa, hat_sigma2) # 修改结果
      } else if (model == "M3") {
        hat_kappa <- mean(E5)
        hat_sigma2 <- ifelse(mean(E6 - 2 * hat_kappa * E5) + hat_kappa^2 > 0, mean(E6 - 2 * hat_kappa * E5) + hat_kappa^2, 0)
        init_para_new[3:4] <- c(hat_kappa, hat_sigma2) # 修改结果
      }
    }
    # 6. gamma ====
    if (model == "M1") {
      nu_hat1 <- tau_hat1 <- list()
      term1 <- as.numeric()
      for (i in 1:n) {
        Lambda_re <- Lambda_fun(type = types, par0 = init_para_new[1:4], t = t_list[[i]], u = u_list[[i]])
        nu_hat1[[i]] <- Lambda_re$Lambda_u
        tau_hat1[[i]] <- Lambda_re$Lambda_t
        term1[i] <- sum(nu_hat1[[i]] + tau_hat1[[i]])
      }
      hat_gamma <- sum(term1) / sum(data, na.rm = T)
      init_para_new[5] <- hat_gamma # 修改结果
    }

    # 终止条件 ====
    if (show_progress) {
      setTxtProgressBar(pb, iter)
    }
    par_re[iter, ] <- init_para_new
    # print(par_re[iter, ])
    # print(sqrt(sum((par_re[iter, ] - par_re[iter-1, ])^2,na.rm=T)))
    if (sqrt(sum((par_re[iter, ] - par_re[iter - 1, ])^2, na.rm = T)) < tol1) break
    # if(max(abs(par_re[iter, ] - par_re[iter-1, ])) < tol1) break
    # print(iter)
  }

  if (model %in% c("M0", "M2")) {
    colnames(par_re) <- c("alpha_t", "beta_t", "alpha_u", "beta_u", "kappa", "sigma2")
  } else if (model == "M1") {
    colnames(par_re) <- c("alpha_t", "beta_t", "alpha_u", "beta_u", "gamma")
  } else if (model == "M3") {
    colnames(par_re) <- c("alpha", "beta", "kappa", "sigma2")
  }
  # 关闭进度条（如果启用了进度条）
  if (show_progress) {
    close(pb)
  }

  return(list(
    "par_re" = par_re,
    "iter" = iter,
    "mse_all" = sqrt(sum((par_re[iter, ] - par_re[iter - 1, ])^2, na.rm = T))
  ))
}


#' Plotting EM Algorithm Iteration Results
#'
#' This function generates a plot showing the parameter estimates over the course of the EM algorithm's iterations.
#' The plot visualizes the convergence of each parameter to its final estimate. If the true values of the parameters
#' are provided, they will be shown as dashed lines on the plot for comparison.
#'
#' The function uses `ggplot2` to create a line plot of the parameter estimates (`par_re`) at each iteration,
#' with an optional comparison to the true parameter values (`par`). The data is reshaped using `pivot_longer`
#' for easier plotting, and the function provides customization options for the plot's appearance.
#'
#' @param em_par A list or data frame containing the EM algorithm's iteration results. Specifically, it should
#'   include the estimated parameters (`par_re`), with rows representing iterations and columns representing
#'   parameters.
#' @param par A vector of true parameter values (optional) to be plotted as horizontal dashed lines for comparison.
#' @param ncol Integer specifying the number of columns in the facet grid layout. The plot will display each parameter
#'   in its own facet.
#' @param ture_value Logical value (`TRUE` or `FALSE`). If `TRUE`, the true parameter values will be plotted as
#'   dashed horizontal lines. If `FALSE`, the true values will not be shown.
#' @param orders A vector of parameter names in the desired order for the facets.
#' @param f_names A vector of labels to be used for each parameter when facetting the plot.
#'
#' @return A `ggplot` object containing the plot of the EM iteration results.
#'
#' @export
#' @import tidyverse ggsci
EM_iter_plot <- function(em_par, par, ncol = 4, ture_value = TRUE,
                         orders, f_names) {
  # em_para$par_re[1:em_para$iter,] EM算法的迭代结果。
  # 添加数学公式
  f_labeller <- function(variable, value) {
    return(f_names[value])
  }

  d1 <- em_par %>%
    data.frame() %>%
    mutate("index" = 1:dim(.)[1]) %>%
    pivot_longer(cols = !index, names_to = "para", values_to = "value")
  d1$para <- factor(d1$para, levels = orders, ordered = TRUE)

  if (ture_value == TRUE) {
    dummy <- data.frame("para" = colnames(em_par), Z = par)
    dummy$para <- factor(dummy$para, levels = orders, ordered = TRUE)
    p1 <- ggplot(d1, aes(index / 100, value, color = para)) +
      geom_line() +
      facet_wrap(vars(para), ncol = ncol, scales = "free", labeller = f_labeller) +
      scale_color_aaas(name = "Parameters") +
      geom_hline(data = dummy, aes(yintercept = Z), linetype = "dashed") +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      xlab("Iteration x 100") +
      ylab("Value")
  } else {
    p1 <- ggplot(d1, aes(index / 100, value, color = para)) +
      geom_line() +
      facet_wrap(vars(para), ncol = ncol, scales = "free", labeller = f_labeller) +
      scale_color_aaas(name = "Parameters") +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "none") +
      xlab("Iteration x 100") +
      ylab("Value")
  }
  # ggsave("figures/EM-plot.pdf", p1, width = 9, height = 6)
  return(p1)
}


#' Compute Log-Likelihood for the M4 Model
#'
#' This function computes the log-likelihood for a specified model (M4 model) using the provided parameter estimates
#' (`hat_para`) and data (`data`). The log-likelihood is calculated based on a statistical model where each data point
#' is assumed to follow a specific distribution defined by the `Lambda_fun` function, which is related to the parameter estimates.
#'
#' The log-likelihood is calculated by iterating through each data column and computing the contribution to the
#' log-likelihood based on the given parameters and the likelihood function for each time series.
#'
#' @param hat_para A vector of parameter estimates. The first two elements correspond to the initial parameters
#'   used in the `Lambda_fun`, and the third element corresponds to another model parameter involved in the
#'   likelihood computation.
#' @param types A vector specifying the types or categories associated with each series, passed to the `Lambda_fun`
#'   function for calculating the time-varying function.
#' @param data A matrix or data frame containing the data points for each series. Rows correspond to observations,
#'   and columns represent different time series.
#' @param t_list A list of time points (or other relevant sequence) for each series, passed to the `Lambda_fun`
#'   for time-specific calculations.
#'
#' @return A numeric value representing the negative log-likelihood of the data under the specified model. The function
#'   returns the negative value to facilitate optimization (minimization).
#'
#' @export
M4.loglik <- function(hat_para, types, data = sim_dat$diff_Y_t, t_list = t_list) {
  m <- dim(data)[1]
  n <- dim(data)[2] # number of data points
  tau_hat <- list()
  for (i in 1:n) {
    Lambda_re <- Lambda_fun(type = types, par0 = hat_para[1:2], t = t_list[[i]], u = NULL)
    tau_hat[[i]] <- Lambda_re$Lambda_t
  }

  term <- matrix(NA, n, m)
  for (i in 1:n) {
    term[i, ] <- -1 / 2 * log(2 * pi) - 3 / 2 * log(data[, i]) + log(tau_hat[[i]]) - 1 / 2 * (hat_para[3] * sqrt(data[, i]) - tau_hat[[i]] / sqrt(data[, i]))^2
  }
  log_like_value <- sum(term, na.rm = T)
  return(-log_like_value)
}


#' Compute Log-Likelihood for Various Models
#'
#' This function calculates the log-likelihood for different models, including "M0", "M1", "M2", "M3", and "M4".
#' The log-likelihood is computed based on the provided data and estimated parameters. The function supports
#' different likelihood formulas for each model type, where the models depend on different parameters such as
#' kappa, sigma2, gamma, and others.
#'
#' @param model A string indicating the model type. It can be one of "M0", "M1", "M2", "M3", or "M4".
#' @param est_par A vector of estimated parameters. These include parameters like kappa, sigma2, and others
#'   specific to each model.
#' @param data A matrix or data frame of size m x n, where m is the number of observations and n is the number of series.
#'   Each column represents a time series.
#' @param t_list A list of time points (or other relevant sequence) for each series, passed to `Lambda_fun`.
#' @param u_list A list of values corresponding to some parameter for each series, passed to `Lambda_fun`.
#' @param types A vector specifying the types or categories for each series, used in the `Lambda_fun` function.
#'
#' @return The computed log-likelihood for the specified model, which is a single numeric value.
#' @export
Log.liklihod <- function(model = model, est_par = em_para_re[[q]], data = coat_dat$diff_Y, t_list = t_list, u_list = u_list, types = types) {
  # data: m x n 矩阵，包含所有的 \tilde{y}_{ij}
  # par: 列表，包含参数 kappa 和 sigma2
  # par1: 列表，包含参数 delta_Lambda 和 delta_Lambda0
  m <- nrow(data)
  n <- ncol(data)
  unit_loglik <- as.numeric()
  Lambda <- list()
  if (model %in% c("M0", "M2")) {
    kappa <- est_par[5]
    if (model == "M0") {
      sigma2 <- est_par[6]
    } else if (model == "M2") {
      sigma2 <- kappa^2 * est_par[6]
    }
    for (i in 1:n) {
      Y_dat <- data[!is.na(data[, i]), i]
      Lambda_re <- Lambda_fun(type = types, par0 = est_par[1:4], t = t_list[[i]], u = u_list[[i]])
      Lambda[[i]] <- Lambda_re$Lambda_t + Lambda_re$Lambda_u
      term1 <- -1 / 2 * log(2 * pi) - 3 / 2 * log(Y_dat) - 1 / 2 * log(1 + Y_dat * sigma2) + log(Lambda[[i]])
      term2 <- -(kappa^2 * Y_dat - 2 * kappa * Lambda[[i]] + (Lambda[[i]]^2) / Y_dat) / (2 * (Y_dat * sigma2 + 1))
      unit_loglik[i] <- sum(term1 + term2)
    }
    loglik <- sum(unit_loglik)
  } else if (model == "M1") {
    gamma <- est_par[5]
    for (i in 1:n) {
      Y_dat <- data[!is.na(data[, i]), i]
      Lambda_re <- Lambda_fun(type = types, par0 = est_par[1:4], t = t_list[[i]], u = u_list[[i]])
      Lambda[[i]] <- Lambda_re$Lambda_t + Lambda_re$Lambda_u
      term1 <- -1 / 2 * log(2 * pi) - 3 / 2 * log(Y_dat) + log(Lambda[[i]])
      term2 <- -1 / 2 * (Lambda[[i]] / sqrt(Y_dat) - gamma * sqrt(Y_dat))^2
      unit_loglik[i] <- sum(term1 + term2)
    }
    loglik <- sum(unit_loglik)
  } else if (model == "M3") {
    kappa <- est_par[3]
    sigma2 <- est_par[4]
    for (i in 1:n) {
      Y_dat <- data[!is.na(data[, i]), i]
      Lambda_re <- Lambda_fun(type = types, par0 = est_par[1:2], t = t_list[[i]], u = NULL)
      Lambda[[i]] <- Lambda_re$Lambda_t
      term1 <- -1 / 2 * log(2 * pi) - 3 / 2 * log(Y_dat) - 1 / 2 * log(1 + Y_dat * sigma2) + log(Lambda[[i]])
      term2 <- -(kappa^2 * Y_dat - 2 * kappa * Lambda[[i]] + (Lambda[[i]]^2) / Y_dat) / (2 * (Y_dat * sigma2 + 1))
      unit_loglik[i] <- sum(term1 + term2)
    }
    loglik <- sum(unit_loglik)
  } else if (model == "M4") {
    gamma <- est_par[3]
    for (i in 1:n) {
      Lambda_re <- Lambda_fun(type = types, par0 = est_par[1:2], t = t_list[[i]], u = NULL)
      Lambda[[i]] <- Lambda_re$Lambda_t
    }
    for (i in 1:n) {
      Y_dat <- data[!is.na(data[, i]), i]
      term1 <- -1 / 2 * log(2 * pi) - 3 / 2 * log(Y_dat) + log(Lambda[[i]])
      term2 <- -1 / 2 * (gamma * sqrt(Y_dat) - Lambda[[i]] / sqrt(Y_dat))^2
      unit_loglik[i] <- sum(term1 + term2)
    }
    loglik <- sum(unit_loglik)
  }
  return(loglik)
}

# optim(rep(1, 6), log_likelihood, type = types,
#       t = t, u = u, data = sim_dat$diff_Y_t,
#       method = "L-BFGS-B",  # 使用 L-BFGS-B 方法以支持参数约束
#       lower = rep(0, 5),    # 设置所有参数的下界为 0
#       upper = rep(Inf, 5),  # 上界可以设置为无穷大
#       control = list(maxit = 3000))$par
