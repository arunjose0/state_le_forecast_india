# 06_holdout_combined_figure.R
# Observed versus predicted e0 for the 2000 and 2010 holdout origins in one figure.

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

in_file <- here("tables", "holdout_validation", "holdout_predictions_main_model_and_benchmarks.csv")
out_dir <- here("figures", "holdout_validation")
if (!file.exists(in_file)) stop("Missing: ", in_file, call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_levels <- c("Bayesian model", "Persistence", "Linear trend")
model_colours <- c("Bayesian model" = "#D55E00", "Persistence" = "#0072B2", "Linear trend" = "#009E73")
model_shapes <- c("Bayesian model" = 16, "Persistence" = 17, "Linear trend" = 15)

scores <- read_csv(in_file, show_col_types = FALSE) |>
  mutate(
    model = recode(
      model,
      "Main model (regularised baseline)" = "Bayesian model",
      "Benchmark: persistence" = "Persistence",
      "Benchmark: linear trend" = "Linear trend"
    )
  ) |>
  filter(model %in% model_levels, is.finite(observed), is.finite(pred_median)) |>
  mutate(
    model = factor(model, levels = model_levels),
    sex = factor(sex, levels = c("Female", "Male")),
    origin = factor(
      paste("Forecast origin", forecast_start_year),
      levels = paste("Forecast origin", sort(unique(forecast_start_year)))
    )
  )

stopifnot(nrow(scores) > 0, n_distinct(scores$forecast_start_year) == 2)

scores |>
  count(origin, sex, model) |>
  print()

axis_limits <- range(c(scores$observed, scores$pred_median)) + c(-1, 1)

p <- ggplot(scores, aes(x = observed, y = pred_median, colour = model, shape = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, colour = "grey45") +
  geom_point(alpha = 0.75, size = 1.5, stroke = 0.25) +
  facet_grid(origin ~ sex) +
  scale_colour_manual(values = model_colours) +
  scale_shape_manual(values = model_shapes) +
  scale_x_continuous(limits = axis_limits) +
  scale_y_continuous(limits = axis_limits) +
  coord_equal() +
  labs(
    x = "Observed life expectancy at birth (years)",
    y = "Predicted life expectancy at birth (years)",
    colour = "Model",
    shape = "Model"
  ) +
  theme_classic(base_size = 10, base_family = "Times New Roman") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(file.path(out_dir, "holdout_observed_vs_predicted_combined.png"), p,
       width = 8, height = 8.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "holdout_observed_vs_predicted_combined.pdf"), p,
       width = 8, height = 8.5, bg = "white", device = cairo_pdf)

cat("Saved: ", file.path(out_dir, "holdout_observed_vs_predicted_combined.png"), "\n", sep = "")
