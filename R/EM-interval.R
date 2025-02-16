#' Bootstrap Confidence Interval for Model Parameters
#'
#' This function applies the bootstrap method to generate a set of simulated samples, estimates model parameters
#' for each sample using the EM algorithm (or MLE), and calculates the confidence intervals for the parameters.
#'
#' @param it Integer, number of bootstrap iterations (default is 100).
#' @param final_par Numeric vector, the final estimated parameters of the model.
#' @param model Character string, model type, one of "M0", "M1", "M2", "M3", or "M4".
#' @param core Integer, number of cores to use for parallel computation.
#' @param types Character string, the type of data (e.g., "type1", "type2", etc.).
#' @param n Number of units.
#' @param m Number of time.
#' @param init_param Numeric vector, initial parameters for the EM algorithm.
#' @param period Numeric vector, the time interval for the analysis, default is `c(0, 1)`.
#' @param t_list Numeric vector, time points.
#' @param u_list Numeric vector, cycle points.
#' @param alpha_value Numeric, significance level for the confidence interval, default is 0.05.
#' @param tol1 Numeric, tolerance for the convergence of the EM algorithm.
#' @param max_iter Integer, maximum number of iterations for the EM algorithm.
#' @param parallel Logical, whether to use parallel computation, default is `TRUE`.
#' @param n_points Integer, the number of points used in the computation.
#' @param n_samples Integer, the number of samples to use in the computation.
#' @param n_intervals Integer, the number of intervals used in the computation.
#' @param show_progress Logical, whether to display the computation progress, default is `FALSE`.
#' @param approx.method Character string, the approximation method to use.
#'
#' @return A list containing two elements:
#'   - `confidence_intervals`: A data frame containing the estimates, standard errors, and lower and upper bounds of the confidence intervals.
#'   - `BT_para`: A matrix containing the parameter estimates from each bootstrap iteration.
#'
#' @export
#' @importFrom stats rnorm na.omit quantile sd
#' @import foreach
#' @import doParallel

CI_Bootstrap <- function(it = 100, final_par, model, core, types, n, m,
                         init_param, period = c(0, 1), t_list, u_list,
                         alpha_value = 0.05, tol1, max_iter, parallel = TRUE,
                         n_points, n_samples, n_intervals,
                         show_progress = FALSE, approx.method) {
  # 初始化矩阵来存储每次迭代的参数
  len1 <- length(final_par)
  BT_para <- matrix(NA, it, len1)
  bootstrap_iter <- function(q) {
    # 为每个自助样本生成随机 gamma (不同模型)
    if (model == "M0") {
      hat_gamma <- rnorm(n, final_par[5], sqrt(final_par[6]))
    } else if (model == "M1") {
      hat_gamma <- final_par[5]
    } else if (model == "M2") {
      hat_gamma <- rnorm(n, 1, sqrt(final_par[6])) # 和M0不同
    } else if (model == "M3") {
      hat_gamma <- rnorm(n, final_par[3], sqrt(final_par[4]))
    } else if (model == "M4") {
      hat_gamma <- final_par[3]
    }
    # 生成自助法模拟数据
    BT_sim_dat <- sim.dat.path(model = model, type = types, m = m, n = n, t = t_list, u = u_list, par = final_par, gamma = hat_gamma)
    # 初始值设定
    inti_par <- init.guess(model = model, types = types, data = BT_sim_dat$diff_Y_t, t_list = t_list, u_list = u_list, init_param = init_param)
    # EM 算法
    if (model %in% c("M0", "M1", "M2", "M3")) {
      em_para <- tryCatch(
        EM(
          data = BT_sim_dat$diff_Y_t, par0 = inti_par, types = types, model = model,
          tol1 = tol1, max_iter = max_iter, t_list = t_list, u_list = u_list, period = period,
          n_points = n_points, n_samples = n_samples, n_intervals = n_intervals,
          approx.method = approx.method, show_progress = show_progress
        ),
        error = function(e) {
          return(NA)
        }
      )

      # 处理 NA 的情况
      if (any(is.na(em_para))) {
        return(rep(NA, len1))
      } else {
        return(as.numeric(em_para$par_re[em_para$iter, ]))
      }
    } else if (model == "M4") {
      # MLE
      em_para_re <- optim(inti_par, M4.loglik,
        types = types, t = t, data = BT_sim_dat$diff_Y_t,
        method = "L-BFGS-B", lower = rep(0, len1), upper = rep(Inf, len1)
      )$par
      return(em_para_re)
    }
  }

  # 并行或非并行执行
  if (parallel) {
    # 设置并行的核心数量
    registerDoParallel(core)
    # 并行循环，使用 %dopar%
    BT_para <- foreach(q = 1:it, .combine = rbind) %dopar% {
      tryCatch(bootstrap_iter(q), error = function(e) {
        return(rep(NA, len1))
      })
    }
    # 停止并行
    stopImplicitCluster()
  } else {
    # 普通循环
    for (q in 1:it) {
      BT_para[q, ] <- tryCatch(bootstrap_iter(q), error = function(e) {
        return(rep(NA, len1))
      })
    }
  }

  # 计算置信区间
  BT_para <- na.omit(BT_para) # dim(BT_para)
  mle_quan_re <- round(apply(BT_para, 2, quantile, c(alpha_value / 2, 0.5, 1 - alpha_value / 2), na.rm = TRUE), 4)

  # 打印结果
  confidence_intervals <- data.frame(
    Estimate = apply(BT_para, 2, mean, na.rm = TRUE),
    StdError = apply(BT_para, 2, sd, na.rm = TRUE),
    Lower = apply(BT_para, 2, quantile, alpha_value / 2, na.rm = TRUE),
    Upper = apply(BT_para, 2, quantile, 1 - alpha_value / 2, na.rm = TRUE)
  )

  return(list(round(confidence_intervals, 3), BT_para))
}
