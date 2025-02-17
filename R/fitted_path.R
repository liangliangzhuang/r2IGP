#' Compute Fitted Path for Each Product
#'
#' This function calculates the fitted degradation values for each product unit based on estimated parameters
#' using the Lambda function. The function assumes two scales, computing the cumulative fitted degradation values
#' for each unit.
#'
#' @param dat A numeric vector of estimated parameters (including Lambda values) for the product.
#' @param select_ID A vector of product IDs to select for the fitted path calculation.
#' @param t_list A list where each element contains the time points for each product.
#' @param u_list A list where each element contains the covariates for each product.
#' @param types A vector specifying the type of model used for the Lambda calculation.
#'
#' @return A tibble containing the fitted degradation values (`Y`) for each product unit over time.
#'
#' @export
fit.path.process <- function(dat, select_ID, t_list, u_list, types) {
  # Initialize empty lists to store the fitted values for each unit
  # 只适合两尺度的
  par_Y <- Yhat_diff_dat <- Yhat_dat <- list()

  for (i in 1:length(t_list)) {
    Lambda_re <- Lambda_fun(type = types, par0 = dat, t = t_list[[i]], u = u_list[[i]])
    Yhat_diff_dat[[i]] <- (Lambda_re$Lambda_t + Lambda_re$Lambda_u) / dat[5]
    Yhat_dat[[i]] <- cumsum(Yhat_diff_dat[[i]])
  }
  Yhat_df <- map2_dfr(Yhat_dat, select_ID, ~ tibble(Unit = .y, Y = .x))
  return(Yhat_df)
}


#' Plot Fitted Path with or without Confidence Intervals
#'
#' This function generates a plot of the fitted degradation paths for each unit, with an option to display confidence
#' intervals based on quantiles. The function supports models "M0", "M1", "M2" and allows for customized plotting.
#'
#' @param dat A numeric matrix containing the estimated parameters for each product.
#' @param True_Y A tibble or data frame containing the true observed values, with columns `TIME`, `UV`, and `observed`.
#' @param t_list A list of time values for each product.
#' @param u_list A list of covariate values for each product (can be NULL).
#' @param select_ID A vector of product IDs to plot.
#' @param types A vector specifying the type of model used for the Lambda calculation.
#' @param ci A logical value indicating whether to include confidence intervals in the plot (default is TRUE).
#' @param scope A numeric vector of length 2 specifying the quantile range for the confidence intervals (e.g., `c(0.025, 0.975)`).
#'
#' @return A ggplot object representing the fitted degradation paths for the selected product units.
#'
#' @export
fit.path.plot <- function(dat, True_Y, t_list, u_list, select_ID, types, ci = TRUE, scope = NULL) {
  # ci 可以根据数据是否包含区间估计，绘制两个版本的结果。
  # 适合 M0，M1,M2
  if (ci == TRUE) {
    fitted_path <- list()
    for (h in 1:dim(dat)[1]) {
      fitted_path[[h]] <- fit.path.process(dat = dat[h, ], t_list = t_list, u_list = u_list, select_ID = select_ID, types = types)
    }
    summmary_fitted_path <- map_dfr(seq_len(nrow(fitted_path[[1]])), function(i) {
      # 提取 fitted_path 中第 i 行所有列表元素的数据
      values <- map_dbl(fitted_path, ~ .x$Y[i])
      # 计算分位数并保留 Unit 信息
      tibble(
        Unit = fitted_path[[1]]$Unit[i],
        Mean = mean(values, na.rm = T),
        Median = median(values, na.rm = T),
        Q_low = quantile(values, probs = scope[1], na.rm = TRUE),
        Q_up = quantile(values, probs = scope[2], na.rm = TRUE)
      )
    })

    Yhat_df <- summmary_fitted_path %>%
      mutate(TIME = True_Y$TIME, UV = True_Y$UV) %>%
      pivot_longer(Mean:Q_up, values_to = "Y", names_to = "Type")
    # 合并
    combined_Y_df <- bind_rows(True_Y, Yhat_df) %>% mutate(Type = factor(Type))
    wide_combined_Y_df <- combined_Y_df %>% pivot_wider(names_from = Type, values_from = Y)
    wide_combined_Y_df_filter <- wide_combined_Y_df
    # 拟合曲线，区间估计！！！
    p1 <- ggplot(data = wide_combined_Y_df_filter, aes(x = TIME)) +
      # 绘制区间估计的阴影区域
      geom_ribbon(aes(ymin = `Q_low`, ymax = `Q_up`), fill = "#D9D9D9", alpha = 0.8) +
      # 绘制观测数据点
      geom_point(aes(y = observed), color = "#3A9B98", size = 0.8) +
      # 绘制均值线
      geom_line(aes(y = Median), color = "#440154", size = 0.8) +
      labs(x = "Time", y = "Degradation") +
      # 使用 facet_wrap 将每个 Unit 单独展示
      # facet_wrap(~ Unit, scales = "free_y",nrow=nrows) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "serif", size = 10),
        axis.text.y = element_text(family = "serif", size = 10),
        legend.title = element_text(family = "serif", size = 12),
        legend.text = element_text(family = "serif", size = 10)
      ) # 设置图例标签的字体)
  } else {
    # 拟合退化值
    Yhat_df <- fit.path.process(dat = dat, t_list = t_list, u_list = u_list, types = types, select_ID = select_ID) %>%
      mutate(TIME = True_Y$TIME, UV = True_Y$UV, Type = "fitted")
    # 合并
    combined_Y_df <- bind_rows(True_Y, Yhat_df) %>% mutate(Type = factor(Type))
    wide_combined_Y_df <- combined_Y_df %>% pivot_wider(names_from = Type, values_from = Y)
    # 不同单元的拟合情况（点估计）
    p1 <- ggplot(wide_combined_Y_df, aes(x = TIME)) +
      # 绘制观测数据点
      geom_point(aes(y = observed), color = "#3A9B98", size = 1.2) +
      # 绘制均值线
      geom_line(aes(y = fitted), color = "#440154", size = 0.8) +
      # 使用 facet_wrap 将每个 Unit 单独展示
      # facet_wrap(vars(Unit),nrow=nrows,scales = "free_y") +
      labs(x = "Time", y = "Degradation") +
      scale_shape_manual(values = c(NA, 1)) + # 设置点的形状和线条
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(family = "serif", size = 10),
        axis.text.y = element_text(family = "serif", size = 10)
      )
  }
  return(p1)
}
