# =============================================================================
# Customisable scatter plot
# =============================================================================
# A reusable plotting helper that produces a scatter plot with optional log
# scaling, density rings, regression lines, and correlation annotations.
#
# Changes from the original create_scatter_plot_v6:
#   - `format_col_name_func` defaults to `format_col_name` (loaded from the
#     pipeline) — no need to pass it in explicitly.
#   - Correlation is computed on the cleaned subset, not the raw input.
#
# Example:
#   source("R/load_all.R")
#   source("analysis/correlations/scatter_plot.R")
#   df <- build_dataset("human")
#   create_scatter_plot(df, x_var = "orfexondensity", y_var = "junctions_density_cds",
#                       color_var = "halflife", add_density_rings = TRUE)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(scales)
})


create_scatter_plot <- function(data,
                                x_var,
                                y_var,
                                color_var,
                                format_col_name_func = format_col_name,
                                log_x = FALSE,
                                log_y = FALSE,
                                log_color = FALSE,
                                add_density_rings = FALSE,
                                quantile_limit = 1,
                                point_alpha = 0.8,
                                point_size = 1,
                                cor_method = c("spearman", "pearson"),
                                point_color_option = "inferno") {
  cor_method <- match.arg(cor_method)

  # --- Input checks -------------------------------------------------------
  if (!is.data.frame(data)) stop("'data' must be a dataframe.")
  missing_vars <- setdiff(c(x_var, y_var, color_var), names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in dataframe: ",
         paste(missing_vars, collapse = ", "))
  }
  if (!is.function(format_col_name_func)) {
    stop("'format_col_name_func' must be a function.")
  }

  # --- Prepare data & labels ---------------------------------------------
  plot_data <- data[!is.na(data[[x_var]]) & !is.na(data[[y_var]]), ]
  dropped <- nrow(data) - nrow(plot_data)
  if (dropped > 0) message("Dropped ", dropped, " rows with missing x/y values.")

  x_lab <- format_col_name_func(x_var)
  y_lab <- format_col_name_func(y_var)
  color_lab <- format_col_name_func(color_var)

  # Auto-discretise numeric colour variables with few unique values
  if (is.numeric(plot_data[[color_var]]) &&
      length(unique(plot_data[[color_var]])) < 10) {
    plot_data[[color_var]] <- as.factor(plot_data[[color_var]])
  }

  # --- Log handling (apply to labels / data, plus axes later) ------------
  if (log_x) {
    if (any(plot_data[[x_var]] <= 0, na.rm = TRUE)) {
      warning("x-variable '", x_var,
              "' contains non-positive values; adding pseudocount of 1.",
              call. = FALSE)
      plot_data[[x_var]] <- plot_data[[x_var]] + 1
    }
    x_lab <- paste0(x_lab, " (log10)")
  }
  if (log_y) {
    if (any(plot_data[[y_var]] <= 0, na.rm = TRUE)) {
      warning("y-variable '", y_var,
              "' contains non-positive values; adding pseudocount of 1.",
              call. = FALSE)
      plot_data[[y_var]] <- plot_data[[y_var]] + 1
    }
    y_lab <- paste0(y_lab, " (log10)")
  }

  # --- Correlation on cleaned data ---------------------------------------
  cor_res <- tryCatch(
    cor.test(plot_data[[x_var]], plot_data[[y_var]],
             method = cor_method, use = "pairwise.complete.obs",
             exact = FALSE),
    error = function(e) NULL
  )
  r_text <- if (!is.null(cor_res)) {
    sprintf("%s R = %.3f, %s",
            tools::toTitleCase(cor_method),
            cor_res$estimate,
            format.pval(cor_res$p.value, digits = 2, eps = 0.001))
  } else {
    "Correlation not available"
  }

  # --- Axis limits -------------------------------------------------------
  x_vals <- plot_data[[x_var]]
  y_vals <- plot_data[[y_var]]
  x_lims <- if (length(x_vals) > 1 && quantile_limit < 1) {
    quantile(x_vals, probs = c((1 - quantile_limit) / 2,
                               1 - (1 - quantile_limit) / 2),
             na.rm = TRUE)
  } else range(x_vals, na.rm = TRUE)
  y_lims <- if (length(y_vals) > 1 && quantile_limit < 1) {
    quantile(y_vals, probs = c((1 - quantile_limit) / 2,
                               1 - (1 - quantile_limit) / 2),
             na.rm = TRUE)
  } else range(y_vals, na.rm = TRUE)

  # --- Plot --------------------------------------------------------------
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(colour = .data[[color_var]]),
               alpha = point_alpha, size = point_size, shape = 16)

  if (add_density_rings) {
    p <- p +
      geom_density_2d(colour = "black", alpha = 0.9, linewidth = 0.6) +
      geom_density_2d(colour = "white", alpha = 0.9, linewidth = 0.3)
  }

  p <- p +
    geom_smooth(method = "lm",    se = FALSE, colour = "#4B0082",
                formula = y ~ x, linewidth = 0.8) +
    geom_smooth(method = "loess", se = FALSE, colour = "#FF4500",
                formula = y ~ x, linewidth = 0.8, span = 0.75) +
    coord_cartesian(xlim = x_lims, ylim = y_lims, expand = TRUE) +
    labs(
      title    = paste(y_lab, "vs.", x_lab),
      subtitle = r_text,
      x = x_lab, y = y_lab, colour = color_lab
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 19, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title    = element_text(size = 18),
      axis.text     = element_text(size = 16),
      legend.position = "right"
    )

  # --- Colour scale ------------------------------------------------------
  is_discrete_colour <- is.factor(plot_data[[color_var]]) ||
                        is.character(plot_data[[color_var]])
  if (is_discrete_colour) {
    p <- p + scale_color_viridis_d(option = point_color_option)
  } else if (log_color) {
    if (any(plot_data[[color_var]] <= 0, na.rm = TRUE)) {
      warning("Colour variable '", color_var,
              "' has non-positive values; shifting for log scale.",
              call. = FALSE)
      shift <- min(plot_data[[color_var]][plot_data[[color_var]] > 0],
                   na.rm = TRUE) / 2
      plot_data[[color_var]] <- plot_data[[color_var]] + shift
    }
    p <- p + scale_color_viridis_c(option = point_color_option,
                                   trans = "log10",
                                   labels = scales::label_log(10))
  } else {
    p <- p + scale_color_viridis_c(option = point_color_option)
  }

  if (log_x) p <- p + scale_x_log10(labels = scales::label_log(10)) +
                      annotation_logticks(sides = "b")
  if (log_y) p <- p + scale_y_log10(labels = scales::label_log(10)) +
                      annotation_logticks(sides = "l")

  p
}
