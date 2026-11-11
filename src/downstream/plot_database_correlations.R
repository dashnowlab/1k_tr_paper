#!/usr/bin/env Rscript

# Compare ATaRVa CN/LPS statistics with TR catalog population databases.
#
# Expected object in the current R session:
#   df_annotated
#
# Source the script in that R/RStudio session:
#
#   source("src/downstream/plot_database_correlations.R")
#
# The script creates:
#   correlation_summary
#   p_lps_correlation
#   p_cn_correlation
#
# It does not display or save plots by default. To save results, change
# save_outputs to TRUE below.

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

maximum_points_per_panel <- 30000L
random_seed <- 20261105L
save_outputs <- FALSE
output_directory <- "database_correlation_plots"

library(dplyr)
library(ggplot2)
library(scales)
library(tibble)

stopifnot(exists("df_annotated"), is.data.frame(df_annotated))

# -----------------------------------------------------------------------------
# Comparison definitions
# -----------------------------------------------------------------------------

# LPS is compared with AoU1027 and HPRC256 LPS values.
# CN is compared with HPRC256 total-length and TenK10K values.
# These comparisons assume that both sides are expressed in repeat-copy units.

comparison_spec <- tibble::tribble(
  ~measurement, ~database, ~metric, ~our_column, ~database_column,
  "LPS", "AoU1027", "StDev", "lps_stdev", "AoU1027_Stdev",
  "LPS", "AoU1027", "Median", "lps_median", "AoU1027_Median",
  "LPS", "AoU1027", "P99", "lps_99th_percentile", "AoU1027_99thPercentile",
  "LPS", "AoU1027", "Maximum", "lps_max", "AoU1027_MaxAllele",
  "LPS", "HPRC256 LPS", "StDev", "lps_stdev", "HPRC256_Stdev",
  "LPS", "HPRC256 LPS", "Median", "lps_median", "HPRC256_Median",
  "LPS", "HPRC256 LPS", "P99", "lps_99th_percentile", "HPRC256_99thPercentile",
  "LPS", "HPRC256 LPS", "Maximum", "lps_max", "HPRC256_MaxAllele",
  "CN", "HPRC256 total length", "StDev", "cn_stdev", "HPRC256TotalLength_Stdev",
  "CN", "HPRC256 total length", "Median", "cn_median", "HPRC256TotalLength_Median",
  "CN", "HPRC256 total length", "P99", "cn_99th_percentile", "HPRC256TotalLength_99thPercentile",
  "CN", "HPRC256 total length", "Maximum", "cn_max", "HPRC256TotalLength_MaxAllele",
  "CN", "TenK10K", "StDev", "cn_stdev", "TenK10K_Stdev",
  "CN", "TenK10K", "Median", "cn_median", "TenK10K_Median",
  "CN", "TenK10K", "P99", "cn_99th_percentile", "TenK10K_99thPercentile",
  "CN", "TenK10K", "Maximum", "cn_max", "TenK10K_MaxAllele"
)

required_columns <- unique(c(
  comparison_spec$our_column,
  comparison_spec$database_column
))
missing_columns <- setdiff(required_columns, names(df_annotated))

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing from df_annotated: ",
    paste(missing_columns, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Numeric conversion helper
# -----------------------------------------------------------------------------

to_numeric <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }

  x <- as.character(x)
  x[x %in% c(".", "", "NA", "NaN", "nan")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

format_correlation <- function(x) {
  ifelse(is.finite(x), sprintf("%.2f", x), "NA")
}

# -----------------------------------------------------------------------------
# Calculate correlations using all matched loci and retain only a plotting
# sample. This avoids constructing a many-million-row table for every metric.
# -----------------------------------------------------------------------------

set.seed(random_seed)

summary_list <- vector("list", nrow(comparison_spec))
plot_list <- vector("list", nrow(comparison_spec))

for (i in seq_len(nrow(comparison_spec))) {
  specification <- comparison_spec[i, ]

  our_values <- to_numeric(
    df_annotated[[specification$our_column]]
  )
  database_values <- to_numeric(
    df_annotated[[specification$database_column]]
  )

  valid_indexes <- which(
    is.finite(our_values) &
      is.finite(database_values) &
      our_values >= 0 &
      database_values >= 0
  )

  matched_loci <- length(valid_indexes)

  if (matched_loci >= 3) {
    pearson_r <- suppressWarnings(
      cor(
        our_values[valid_indexes],
        database_values[valid_indexes],
        method = "pearson"
      )
    )

    spearman_rho <- suppressWarnings(
      cor(
        our_values[valid_indexes],
        database_values[valid_indexes],
        method = "spearman"
      )
    )
  } else {
    pearson_r <- NA_real_
    spearman_rho <- NA_real_
  }

  summary_list[[i]] <- data.frame(
    measurement = specification$measurement,
    database = specification$database,
    metric = specification$metric,
    our_column = specification$our_column,
    database_column = specification$database_column,
    matched_loci = matched_loci,
    pearson_r = pearson_r,
    spearman_rho = spearman_rho,
    stringsAsFactors = FALSE
  )

  if (matched_loci > maximum_points_per_panel) {
    plot_indexes <- sample(valid_indexes, maximum_points_per_panel)
  } else {
    plot_indexes <- valid_indexes
  }

  plotted_loci <- length(plot_indexes)

  plot_list[[i]] <- data.frame(
    measurement = rep(specification$measurement, plotted_loci),
    database = rep(specification$database, plotted_loci),
    metric = rep(specification$metric, plotted_loci),
    our_value = our_values[plot_indexes],
    database_value = database_values[plot_indexes],
    stringsAsFactors = FALSE
  )

  rm(our_values, database_values, valid_indexes, plot_indexes)

  if (i %% 4 == 0) {
    invisible(gc())
  }
}

correlation_summary <- bind_rows(summary_list) %>%
  mutate(
    label = paste0(
      "n = ", scales::comma(matched_loci),
      "\nPearson r = ", format_correlation(pearson_r),
      "\nSpearman ρ = ", format_correlation(spearman_rho)
    )
  )

plot_data <- bind_rows(plot_list)

print(
  correlation_summary %>%
    select(
      measurement,
      database,
      metric,
      matched_loci,
      pearson_r,
      spearman_rho
    )
)

# -----------------------------------------------------------------------------
# Plot function
# -----------------------------------------------------------------------------

make_correlation_plot <- function(measurement_name) {
  points <- plot_data %>%
    filter(measurement == measurement_name)

  annotations <- correlation_summary %>%
    filter(measurement == measurement_name)

  ggplot(
    points,
    aes(
      x = our_value,
      y = database_value
    )
  ) +
    geom_point(
      alpha = 0.15,
      size = 0.45,
      color = "#0072B2"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "#D55E00",
      linewidth = 0.7
    ) +
    geom_label(
      data = annotations,
      aes(
        x = -Inf,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.05,
      size = 3,
      label.size = 0.2,
      fill = scales::alpha("white", 0.88)
    ) +
    facet_grid(
      rows = vars(metric),
      cols = vars(database),
      scales = "free"
    ) +
    scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1)
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1)
    ) +
    labs(
      title = paste0(
        measurement_name,
        " comparison with external databases"
      ),
      subtitle = paste0(
        "Correlations use all matched loci; plots show up to ",
        scales::comma(maximum_points_per_panel),
        " loci per panel"
      ),
      x = paste0("Our ", measurement_name, " value"),
      y = "External database value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92"),
      plot.title.position = "plot"
    )
}

p_lps_correlation <- make_correlation_plot("LPS")
p_cn_correlation <- make_correlation_plot("CN")

# Display manually when ready:
# print(p_lps_correlation)
# print(p_cn_correlation)

# -----------------------------------------------------------------------------
# Optional output
# -----------------------------------------------------------------------------

if (isTRUE(save_outputs)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  write.csv(
    correlation_summary,
    file.path(output_directory, "database_correlation_summary.csv"),
    row.names = FALSE
  )

  ggsave(
    filename = file.path(output_directory, "lps_database_correlations.png"),
    plot = p_lps_correlation,
    width = 11,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(output_directory, "cn_database_correlations.png"),
    plot = p_cn_correlation,
    width = 11,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  message(
    "Correlation results written to: ",
    normalizePath(output_directory)
  )
}
