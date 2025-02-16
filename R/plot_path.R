#' Degradation Path Plot Summary
#' This function generates the degradation path plot based on different lambda functions for time and cycles.
#'
#' @param data List of degradation data.
#' @param model Character string specifying the model type (e.g., "M0", "M1", etc.).
#' @param t Numeric vector of time.
#' @param u Numeric vector for units or cycles.
#' @param leg.pos Character string for legend position. Default is "none".
#'
#' @return A list of ggplot objects.
#' @import ggplot2 dplyr forcats viridis
#' @export
degradation.path.plot.summary <- function(data = sim_dat$diff_Y_t, model, t = t, u = NULL, leg.pos = "none") {
  n <- dim(data[[1]])[2]

  # Prepare data by adding a row of zeros to simulate starting point
  sim_dat_plot <- lapply(data, function(df) {
    new_row <- rep(0, n)
    rbind(new_row, df)
  })

  # For models "M0", "M1", and "M2"
  if (model %in% c("M0", "M1", "M2")) {
    # Transform the data to data frames and label groups (X, Y, Z)
    df1 <- as.data.frame(sim_dat_plot[[1]]) %>% mutate(group = "Y", time = t)
    df2 <- as.data.frame(sim_dat_plot[[2]]) %>% mutate(group = "X", time = t)
    df3 <- as.data.frame(sim_dat_plot[[3]]) %>% mutate(group = "Z", time = u)

    # Combine X and Z groups for plotting
    df_all <- bind_rows(df2, df3)
    new_df_all <- df_all %>% pivot_longer(cols = starts_with("n"), names_to = "n", values_to = "value")

    # Adjust for secondary x-axis scaling
    ind1 <- u[length(u)] / t[length(t)]
    t_breaks <- seq(min(t), max(t), length.out = 5)
    u_breaks <- round(u[t_breaks + 1], 2)

    # Plot the X(t) and Z(t) paths
    p1 <- ggplot() +
      geom_line(data = new_df_all %>% filter(group == "X"), aes(x = time, y = value, linetype = n), color = "#64A4A3", alpha = 1) +
      geom_point(data = new_df_all %>% filter(group == "X"), aes(x = time, y = value), color = "#64A4A3", size = 1) +
      geom_line(data = new_df_all %>% filter(group == "Z"), aes(x = time / ind1, y = value, linetype = n), color = "#5D336F", alpha = 1) +
      geom_point(data = new_df_all %>% filter(group == "Z"), aes(x = time / ind1, y = value), color = "#5D336F", size = 1) +
      scale_x_continuous(name = "t", sec.axis = sec_axis(~ . * ind1, name = "u")) +
      theme_bw() +
      ylab("X(t) or Z(t)") +
      scale_linetype_manual(name = "", values = rep(1, n)) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = leg.pos,
        axis.line.x.bottom = element_line(color = "#64A4A3"),
        axis.line.x.top = element_line(color = "#5D336F"),
        axis.text.x.bottom = element_text(color = "#64A4A3", family = "serif", size = 10),
        axis.text.x.top = element_text(color = "#5D336F", family = "serif", size = 10),
        axis.title.x.bottom = element_text(color = "#64A4A3", family = "serif", size = 12),
        axis.title.x.top = element_text(color = "#5D336F", family = "serif", size = 12),
        axis.text.y = element_text(family = "serif", size = 10),
        axis.title.y = element_text(family = "serif", size = 12),
        legend.title = element_text(family = "serif", size = 12)
      )

    # Plot the Y(t) paths
    df_Y <- df1 %>%
      pivot_longer(cols = starts_with("n"), names_to = "n", values_to = "value") %>%
      mutate(t = time)
    p2 <- ggplot(df_Y, aes(x = t, y = value, color = forcats::fct(n))) +
      geom_line() +
      geom_point() +
      scale_x_continuous(name = "t", sec.axis = sec_axis(~ . * ind1, name = "u")) +
      scale_color_viridis(discrete = TRUE) +
      ylab("Y(t)") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = leg.pos,
        axis.line.x.bottom = element_line(color = "#64A4A3"),
        axis.line.x.top = element_line(color = "#5D336F"),
        axis.text.x.bottom = element_text(color = "#64A4A3", family = "serif", size = 10),
        axis.text.x.top = element_text(color = "#5D336F", family = "serif", size = 10),
        axis.title.x.bottom = element_text(color = "#64A4A3", family = "serif", size = 12),
        axis.title.x.top = element_text(color = "#5D336F", family = "serif", size = 12),
        axis.text.y = element_text(family = "serif", size = 10),
        axis.title.y = element_text(family = "serif", size = 12),
        legend.title = element_text(family = "serif", size = 12)
      )

    return(list(p1, p2))
  } else if (model %in% c("M3", "M4")) {
    df1 <- as.data.frame(sim_dat_plot[[1]]) %>% mutate(group = "Y", t = t)
    df_Y <- df1 %>% pivot_longer(!c(group, t), values_to = "value", names_to = "n")
    p2 <- ggplot(df_Y, aes(t, value, color = forcats::fct(n))) +
      geom_line() +
      geom_point() +
      scale_color_viridis(discrete = TRUE) +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = leg.pos)
    return(p2)
  }
}


#' 3D Path Plot of Degradation Loss
#'
#' Creates a 3D plot to visualize degradation loss for different units across time and cycles.
#'
#' @param model Character string specifying the model type (e.g., "M0", "M1", "M2").
#' @param units Integer specifying the number of units to plot.
#' @param time Numeric vector representing the time data.
#' @param cycles Numeric vector representing the cycle data.
#' @param loss Matrix of loss values, typically from `sim_dat$Y_t`.
#'
#' @return A 3D plot using `plot_ly`.
#' @import plotly RColorBrewer
#' @export
path.3D.plot <- function(model, units = 10, time = t, cycles = u, loss = sim_dat$Y_t) {
  if (model %in% c("M0", "M1", "M2")) {
    colors <- brewer.pal(10, "Set3")
    n <- dim(loss)[2]
    loss <- rbind(rep(0, n), loss) # Add zero to start the plot

    # Initialize the 3D plot
    p <- plot_ly()

    # Add trace for each unit
    for (i in 1:units) {
      p <- add_trace(p,
        x = cycles,
        y = time,
        z = loss[, i],
        type = "scatter3d",
        mode = "lines",
        name = paste("Unit", i),
        line = list(width = 2, color = colors[i])
      )

      # Projection to x-z plane
      p <- add_trace(p,
        x = cycles,
        y = rep(0, length(time)),
        z = loss[, i],
        type = "scatter3d",
        mode = "lines",
        showlegend = FALSE,
        line = list(color = "grey", dash = "dash")
      )

      # Projection to y-z plane
      p <- add_trace(p,
        x = rep(0, length(cycles)),
        y = time,
        z = loss[, i],
        type = "scatter3d",
        mode = "lines",
        showlegend = FALSE,
        line = list(color = "grey", dash = "dash")
      )

      # Projection to x-y plane
      p <- add_trace(p,
        x = cycles,
        y = time,
        z = rep(0, length(loss[, i])),
        type = "scatter3d",
        mode = "lines",
        showlegend = FALSE,
        line = list(color = "black", dash = "dash")
      )
    }

    # Customize layout with axis labels and camera angle
    p <- layout(p,
      scene = list(
        xaxis = list(title = "Cycles", range = c(0, max(cycles))),
        yaxis = list(title = "Time (Days)", range = c(0, max(time))),
        zaxis = list(title = "Capacity Loss (x 100%)"),
        camera = list(eye = list(x = -1, y = -1.5, z = 2))
      ),
      legend = list(x = 0.7, y = 0.7, xanchor = "center", yanchor = "middle")
    )

    return(p)
  } else {
    return("The input data dimensions are incorrect.")
  }
}
