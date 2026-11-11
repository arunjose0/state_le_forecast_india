# Compare the three prior settings using only two figures:
# 1) observed history and forecast trajectories for both sexes in one state-faceted graph;
# 2) posterior-predictive intervals for all three prior settings, faceted by state.
# Run after the main, tighter-prior, and wide-prior model fits have passed diagnostics.

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(rstan)
})

model_dirs <- c(
  "2" = "main_regularised",
  "3" = "tighter_prior_sensitivity",
  "4" = "wide_prior_sensitivity"
)
labels <- c(`2` = "Main model (regularised baseline)",
            `3` = "Tighter-prior sensitivity",
            `4` = "Wide-prior sensitivity")
model_colours <- c(
  "Main model (regularised baseline)" = "#D55E00",
  "Tighter-prior sensitivity" = "#0072B2",
  "Wide-prior sensitivity" = "#009E73"
)
fit_paths <- c(
  "2" = here::here("fits", "fit_main_regularised.rds"),
  "3" = here::here("fits", "fit_tighter_priors.rds"),
  "4" = here::here("fits", "fit_wide_priors.rds")
)
input_paths <- file.path(here::here("tables"), unname(model_dirs), "forecast_summary.csv")
if (any(!file.exists(input_paths))) {
  stop("Missing forecast summary table(s):\n", paste(input_paths[!file.exists(input_paths)], collapse = "\n"),
       "\nRun codes/01_run_models.R and confirm all posterior fits passed diagnostics.", call. = FALSE)
}
if (any(!file.exists(fit_paths))) {
  stop("Missing saved posterior fit(s) needed for the comparative PPC:\n",
       paste(fit_paths[!file.exists(fit_paths)], collapse = "\n"), call. = FALSE)
}
prep_path <- here::here("fits", "prep_objects.rds")
if (!file.exists(prep_path)) stop("Missing preparation object: ", prep_path, call. = FALSE)

output_tables <- here::here("tables", "prior_sensitivity_comparison")
output_figures <- here::here("figures", "prior_sensitivity_comparison")
dir.create(output_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figures, recursive = TRUE, showWarnings = FALSE)

# Remove obsolete figures from earlier runs so this folder contains only the two
# current sensitivity figures after the script completes.
expected_figures <- c(
  "forecast_trajectories_both_sexes_prior_sensitivity.png",
  "posterior_predictive_check_by_state_prior_sensitivity.png"
)
old_pngs <- list.files(output_figures, pattern = "\\.png$", full.names = TRUE)
stale_pngs <- old_pngs[!basename(old_pngs) %in% expected_figures]
if (length(stale_pngs)) unlink(stale_pngs)

summary_all <- purrr::map2_dfr(input_paths, names(model_dirs), function(path, pm) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(prior_mode = as.integer(pm), model_label = unname(labels[pm]))
}) %>%
  mutate(
    forecast_type = factor(forecast_type, levels = c("Latent trajectory", "Future observation")),
    sex = factor(sex, levels = c("Female", "Male")),
    model_label = factor(model_label, levels = unname(labels))
  )
write_csv(summary_all, file.path(output_tables, "forecast_summaries_all_prior_settings.csv"))

# Compact state-year table for transparent numerical sensitivity assessment.
compare <- summary_all %>%
  select(forecast_type, sex, state, g, year, prior_mode,
         median, q025, q05, q25, q75, q95, q975, width50, width90, width95) %>%
  pivot_wider(
    names_from = prior_mode,
    values_from = c(median, q025, q05, q25, q75, q95, q975, width50, width90, width95),
    names_glue = "{.value}_mode{prior_mode}"
  ) %>%
  mutate(
    shift_tight_vs_main = median_mode3 - median_mode2,
    shift_wide_vs_main = median_mode4 - median_mode2,
    width95_change_tight = width95_mode3 - width95_mode2,
    width95_change_wide = width95_mode4 - width95_mode2,
    interval95_overlap_tight = pmax(q025_mode2, q025_mode3) <= pmin(q975_mode2, q975_mode3),
    interval95_overlap_wide = pmax(q025_mode2, q025_mode4) <= pmin(q975_mode2, q975_mode4)
  ) %>%
  arrange(forecast_type, sex, year, state)
write_csv(compare, file.path(output_tables, "state_year_numeric_sensitivity_comparison.csv"))

horizons <- intersect(c(2030, 2040, 2050), unique(summary_all$year))
width_summary <- summary_all %>% filter(year %in% horizons) %>%
  group_by(forecast_type, sex, year, model_label) %>%
  summarise(
    n_states = n(),
    median_width50 = median(width50, na.rm = TRUE),
    median_width90 = median(width90, na.rm = TRUE),
    median_width95 = median(width95, na.rm = TRUE),
    q25_width95 = quantile(width95, .25, na.rm = TRUE),
    q75_width95 = quantile(width95, .75, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(width_summary, file.path(output_tables, "summary_interval_widths_by_horizon.csv"))

shift_summary <- compare %>% filter(year %in% horizons) %>%
  select(forecast_type, sex, state, year, shift_tight_vs_main, shift_wide_vs_main) %>%
  pivot_longer(c(shift_tight_vs_main, shift_wide_vs_main),
               names_to = "comparison", values_to = "median_shift_years") %>%
  mutate(comparison = recode(comparison,
    shift_tight_vs_main = "Tighter prior minus main",
    shift_wide_vs_main = "Wider prior minus main")) %>%
  group_by(forecast_type, sex, year, comparison) %>%
  summarise(
    n_states = n(),
    median_absolute_shift = median(abs(median_shift_years), na.rm = TRUE),
    mean_absolute_shift = mean(abs(median_shift_years), na.rm = TRUE),
    maximum_absolute_shift = max(abs(median_shift_years), na.rm = TRUE),
    states_shift_over_1_year = sum(abs(median_shift_years) > 1, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(shift_summary, file.path(output_tables, "summary_median_shifts_by_horizon.csv"))

prep <- readRDS(prep_path)
df <- tibble::as_tibble(prep$df)
stan_base <- prep$stan_data_full
if (nrow(df) != stan_base$N) stop("prep_objects.rds rows do not align with Stan data.", call. = FALSE)

# Figure 1: one combined-sex trajectory graph, faceted by state.
history <- df %>%
  transmute(state = as.character(state), year = as.numeric(year),
            Female = as.numeric(e0_f), Male = as.numeric(e0_m)) %>%
  pivot_longer(c(Female, Male), names_to = "sex", values_to = "observed") %>%
  filter(is.finite(observed))
lat <- summary_all %>% filter(forecast_type == "Latent trajectory") %>%
  mutate(state = as.character(state), sex = as.character(sex), model_label = as.character(model_label))
state_levels <- unique(as.character(prep$state_levels_full))
history <- history %>% mutate(state = factor(state, levels = state_levels))
lat <- lat %>% mutate(state = factor(state, levels = state_levels))

p_trajectory <- ggplot() +
  geom_ribbon(data = lat,
              aes(x = year, ymin = q05, ymax = q95, fill = model_label,
                  group = interaction(model_label, sex)),
              alpha = .08, colour = NA) +
  geom_line(data = lat,
            aes(x = year, y = median, colour = model_label, linetype = sex,
                group = interaction(model_label, sex)),
            linewidth = .55) +
  geom_point(data = history,
             aes(x = year, y = observed, shape = sex),
             inherit.aes = FALSE, colour = "grey35", size = .65, alpha = .7) +
  facet_wrap(~state, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = model_colours) +
  scale_fill_manual(values = model_colours) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  labs(x = "Year", y = "Life expectancy at birth (years)",
       colour = "Prior setting", fill = "Prior setting", linetype = "Sex", shape = "Observed sex") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(size = 8), panel.spacing = unit(.5, "lines"))
ggsave(file.path(output_figures, "forecast_trajectories_both_sexes_prior_sensitivity.png"),
       p_trajectory, width = 16, height = 17, dpi = 600, bg = "white")

# Figure 2: one comparative posterior-predictive interval plot, faceted by state.
# The first observation per state fixes the model's starting level and is excluded,
# because its generated y_rep value is intentionally set to the observed value.
include <- stan_base$is_first == 0L
idx <- which(include)
if (!length(idx)) stop("No non-initial observations available for comparative PPC.", call. = FALSE)
ppc_rows <- list()
ppc_coverage_rows <- list()
row_k <- 0L
for (pm in as.integer(names(fit_paths))) {
  model_label <- unname(labels[as.character(pm)])
  cat("Extracting posterior-predictive draws for ", model_label, "\n", sep = "")
  fit <- readRDS(unname(fit_paths[as.character(pm)]))
  ex <- rstan::extract(fit, pars = c("y_rep_f", "y_rep_m"), permuted = TRUE)
  for (sx in c("Female", "Male")) {
    observed <- if (sx == "Female") as.numeric(df$e0_f[idx]) else as.numeric(df$e0_m[idx])
    reps <- if (sx == "Female") ex$y_rep_f[, idx, drop = FALSE] else ex$y_rep_m[, idx, drop = FALSE]
    if (!all(is.finite(reps))) {
      stop("Non-finite posterior-predictive draws for ", model_label, " / ", sx, call. = FALSE)
    }
    qs <- t(apply(reps, 2, quantile, probs = c(.05, .25, .75, .95),
                  names = FALSE, type = 8, na.rm = FALSE))
    row_k <- row_k + 1L
    ppc_rows[[row_k]] <- tibble(
      state = as.character(df$state[idx]), year = as.numeric(df$year[idx]), sex = sx,
      observed = observed, q05 = qs[, 1], q25 = qs[, 2], q75 = qs[, 3], q95 = qs[, 4],
      model_label = model_label, prior_mode = pm
    )
    ppc_coverage_rows[[row_k]] <- tibble(
      model_label = model_label, prior_mode = pm, sex = sx,
      n_observations = length(observed),
      coverage50 = mean(observed >= qs[, 2] & observed <= qs[, 3]),
      coverage90 = mean(observed >= qs[, 1] & observed <= qs[, 4])
    )
  }
  rm(ex, fit)
  invisible(gc(FALSE))
}
ppc_data <- bind_rows(ppc_rows) %>%
  mutate(state = factor(state, levels = state_levels),
         model_label = factor(model_label, levels = unname(labels)),
         sex = factor(sex, levels = c("Female", "Male")))
observed_data <- ppc_data %>% distinct(state, year, sex, observed)
position_ppc <- position_dodge(width = .65)
p_ppc <- ggplot(ppc_data,
                aes(x = year, colour = model_label,
                    group = interaction(model_label, sex))) +
  geom_linerange(aes(ymin = q05, ymax = q95, linetype = sex),
                 position = position_ppc, linewidth = .35, alpha = .6) +
  geom_linerange(aes(ymin = q25, ymax = q75, linetype = sex),
                 position = position_ppc, linewidth = .8) +
  geom_point(data = observed_data,
             aes(x = year, y = observed, shape = sex),
             inherit.aes = FALSE, colour = "black", size = .65, alpha = .8) +
  facet_wrap(~state, ncol = 4, scales = "free_x") +
  scale_colour_manual(values = model_colours) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  labs(x = "Year", y = "Life expectancy at birth (years)",
       colour = "Prior setting", linetype = "Sex", shape = "Observed sex") +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(size = 8), panel.spacing = unit(.5, "lines"))
ggsave(file.path(output_figures, "posterior_predictive_check_by_state_prior_sensitivity.png"),
       p_ppc, width = 16, height = 17, dpi = 600, bg = "white")
write_csv(bind_rows(ppc_coverage_rows),
          file.path(output_tables, "posterior_predictive_coverage_by_prior_and_sex.csv"))

cat("\nPRIOR SENSITIVITY COMPARISON COMPLETE\n")
cat("Exactly two comparison figures saved under: ", output_figures, "\n", sep = "")
cat("  1. forecast_trajectories_both_sexes_prior_sensitivity.png\n")
cat("  2. posterior_predictive_check_by_state_prior_sensitivity.png\n")
cat("Numerical summaries saved under: ", output_tables, "\n", sep = "")
