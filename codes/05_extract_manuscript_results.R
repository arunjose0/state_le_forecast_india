# 05_extract_manuscript_results.R
#
# Extract all numerical results needed for the planned manuscript Results section.
# This script DOES NOT fit or refit any model. It reads the prepared data, saved
# posterior fits, forecast tables, holdout results, threshold-probability output,
# diagnostics, posterior predictive checks, and prior-sensitivity tables.
#
# Results sections covered:
#   1. Historical interstate patterns (1971 and 2024; Health Index groups)
#   2. State-level projections to 2050, interstate dispersion, sex gaps,
#      state-specific equilibrium targets, and P(e0 < 80) in 2024 vs 2050
#   3. Model-estimated COVID deviations and post-2021 decay half-lives
#   4. Out-of-sample validation metrics (origins 2000 and 2010)
#   5. Sampling diagnostics, posterior/prior predictive checks, and prior sensitivity
#
# Run from the project root, after scripts 00-04 and prior-sensitivity analysis:
#   Rscript codes/05_extract_manuscript_results.R
#
# Main output folder:
#   tables/manuscript_results_raw/

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
})

# -----------------------------------------------------------------------------
# 1. Project paths and inputs
# -----------------------------------------------------------------------------

project_root <- here::here()
out_dir <- file.path(project_root, "tables", "manuscript_results_raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Required inputs for the substantive Results sections.
required_rel <- c(
  "fits/prep_objects.rds",
  "fits/fit_main_regularised.rds",
  "tables/main_regularised/forecast_summary.csv",
  "tables/main_regularised/state_equilibrium_targets.csv",
  "tables/main_regularised/mcmc_summary_all_quantities.csv",
  "tables/main_regularised/diagnostic_status.csv",
  "tables/main_regularised/generated_quantities_status.csv",
  "tables/main_regularised/historical_one_step_ppc_coverage.csv",
  "tables/holdout_validation/holdout_metrics_overall.csv",
  "tables/holdout_validation/holdout_metrics_by_origin_sex_lead.csv",
  "tables/threshold_probabilities_2024_2050/state_probability_below_80_2024_2050.csv",
  "tables/threshold_probabilities_2024_2050/health_index_group_probability_below_80_2024_2050.csv"
)

# Optional inputs enrich diagnostics and sensitivity reporting; their absence is
# reported in an inventory and does not silently imply a successful result.
optional_rel <- c(
  "tables/main_regularised/sampler_diagnostics.csv",
  "tables/main_regularised/bfmi_by_chain.csv",
  "tables/main_regularised/ppc_observation_index_key.csv",
  "tables/holdout_validation/holdout_metrics_by_state.csv",
  "tables/holdout_validation/holdout_predictions_main_model_and_benchmarks.csv",
  "tables/holdout_validation/sampler_diagnostics_main_model_holdout_fits.csv",
  "tables/holdout_validation/states_without_pre_origin_training_data.csv",
  "tables/prior_predictive/prior_predictive_historical_interval_coverage.csv",
  "tables/prior_predictive/prior_predictive_plausibility_shares.csv",
  "tables/prior_predictive/prior_predictive_quantile_ranges_001_099.csv",
  "tables/prior_predictive/prior_predictive_forecast_quantiles_all_states_years.csv",
  "tables/prior_sensitivity_comparison/state_year_numeric_sensitivity_comparison.csv",
  "tables/prior_sensitivity_comparison/summary_interval_widths_by_horizon.csv",
  "tables/prior_sensitivity_comparison/summary_median_shifts_by_horizon.csv",
  "tables/prior_sensitivity_comparison/posterior_predictive_coverage_by_prior_and_sex.csv",
  "tables/main_regularised/forecast_summary_last_year.csv",
  "tables/main_regularised/key_parameter_summaries.csv"
)

model_dirs <- c(
  main_regularised = "Main model (regularised baseline)",
  tighter_prior_sensitivity = "Tighter-prior sensitivity",
  wide_prior_sensitivity = "Wide-prior sensitivity"
)

# Required files for sensitivity/diagnostic comparison, but treated as optional
# as a group so the extraction can still produce the core Results if a sensitivity
# run has not yet completed. Missing files are recorded below.
for (d in names(model_dirs)) {
  optional_rel <- c(optional_rel,
    file.path("tables", d, "diagnostic_status.csv"),
    file.path("tables", d, "generated_quantities_status.csv"),
    file.path("tables", d, "mcmc_summary_all_quantities.csv"),
    file.path("tables", d, "forecast_summary.csv"),
    file.path("tables", d, "historical_one_step_ppc_coverage.csv")
  )
}
for (pm in 2:4) {
  optional_rel <- c(optional_rel,
    file.path("tables", "prior_predictive", paste0("mode_", pm), "diagnostic_status.csv"),
    file.path("tables", "prior_predictive", paste0("mode_", pm), "generated_quantities_status.csv")
  )
}
optional_rel <- unique(optional_rel)

input_inventory <- tibble(
  path = unique(c(required_rel, optional_rel)),
  required = path %in% required_rel
) %>%
  mutate(exists = file.exists(file.path(project_root, path)))
write_csv(input_inventory, file.path(out_dir, "input_file_inventory.csv"))

missing_required <- input_inventory %>% filter(required, !exists) %>% pull(path)
if (length(missing_required)) {
  stop(
    "Required input(s) missing under project root: ", project_root,
    "\n - ", paste(missing_required, collapse = "\n - "),
    "\nRun the data preparation, main model, holdout validation, and threshold-probability script first.",
    call. = FALSE
  )
}
missing_optional <- input_inventory %>% filter(!required, !exists) %>% pull(path)
if (length(missing_optional)) {
  message("Some optional diagnostic/sensitivity files are not present. They will be marked unavailable in input_file_inventory.csv.")
}

read_required <- function(rel, name_repair = "unique") {
  path <- file.path(project_root, rel)
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, name_repair = name_repair)
}
read_optional <- function(rel, name_repair = "unique") {
  path <- file.path(project_root, rel)
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, name_repair = name_repair)
}
write_table <- function(x, filename) {
  readr::write_csv(x, file.path(out_dir, filename), na = "")
  message("Wrote: ", filename)
  invisible(x)
}

# -----------------------------------------------------------------------------
# 2. Read source data and apply the same Health Index group labels as the model
# -----------------------------------------------------------------------------

prep <- readRDS(file.path(project_root, "fits", "prep_objects.rds"))
df <- tibble::as_tibble(prep$df) %>%
  transmute(
    state = as.character(state),
    year = as.numeric(year),
    g = as.integer(g),
    e0_f = as.numeric(e0_f),
    e0_m = as.numeric(e0_m)
  )
# Preserve prep$df row order: Stan data, is_first flags, and generated y_rep draws
# use that exact order. Ordering is applied locally only where needed below.
state_names <- as.character(prep$state_levels_full)
stan_data <- prep$stan_data_full
main_fit <- readRDS(file.path(project_root, "fits", "fit_main_regularised.rds"))

if (nrow(df) != stan_data$N) stop("prep$df and Stan data N do not align.", call. = FALSE)
if (length(state_names) != stan_data$N_groups) stop("State labels do not align with Stan groups.", call. = FALSE)
if (!all(c("state", "year", "g", "e0_f", "e0_m") %in% names(df))) {
  stop("Prepared data must contain state, year, g, e0_f, and e0_m.", call. = FALSE)
}

state_tier_map <- tribble(
  ~state_name,              ~health_index_group,
  "Andhra Pradesh",         "Front-runner",
  "Telangana",              "Front-runner",
  "Kerala",                 "Front-runner",
  "Maharashtra",            "Front-runner",
  "Tamil Nadu",             "Front-runner",
  "West Bengal",            "Front-runner",
  "Chhattisgarh",           "Achiever",
  "Goa",                    "Achiever",
  "Gujarat",                "Achiever",
  "Haryana",                "Achiever",
  "Himachal Pradesh",       "Achiever",
  "Karnataka",              "Achiever",
  "Punjab",                 "Achiever",
  "Assam",                  "Aspirant",
  "Bihar",                  "Aspirant",
  "Delhi",                  "Aspirant",
  "Jammu & Kashmir",        "Aspirant",
  "Jharkhand",              "Aspirant",
  "Madhya Pradesh",         "Aspirant",
  "Odisha",                 "Aspirant",
  "Rajasthan",              "Aspirant",
  "Uttar Pradesh",          "Aspirant",
  "Uttarakhand",            "Aspirant"
)

normalise_state <- function(x) {
  z <- tolower(gsub("[^[:alnum:]]", "", trimws(as.character(x))))
  dplyr::case_when(
    z == "chhatisgarh" ~ "chhattisgarh",
    z == "orissa" ~ "odisha",
    TRUE ~ z
  )
}
state_tier_map <- state_tier_map %>%
  mutate(state_key = normalise_state(state_name)) %>%
  distinct(state_key, .keep_all = TRUE)
state_group_index <- match(normalise_state(state_names), state_tier_map$state_key)
if (anyNA(state_group_index)) {
  stop("No Health Index group mapped for: ", paste(state_names[is.na(state_group_index)], collapse = ", "), call. = FALSE)
}
state_groups <- state_tier_map$health_index_group[state_group_index]
state_group_lookup <- tibble(
  state = state_names,
  health_index_group = state_groups,
  g = seq_along(state_names)
)
write_table(state_group_lookup, "state_health_index_group_key.csv")

history_long <- df %>%
  pivot_longer(c(e0_f, e0_m), names_to = "sex_variable", values_to = "e0") %>%
  mutate(sex = recode(sex_variable, e0_f = "Female", e0_m = "Male")) %>%
  left_join(state_group_lookup %>% select(state, health_index_group), by = "state") %>%
  filter(is.finite(year), is.finite(e0))

# State-by-sex values for the two historical endpoints, used to draft ranked
# top-five and bottom-five summaries directly from the extracted-results folder.
historical_state_values <- history_long %>%
  filter(year %in% c(1971, 2024)) %>%
  select(year, sex, state, health_index_group, e0) %>%
  arrange(year, sex, e0, state)
write_table(historical_state_values, "historical_state_values_1971_2024.csv")

# -----------------------------------------------------------------------------
# 3. Results Section 1: historical interstate patterns
# -----------------------------------------------------------------------------

hist_years <- c(1971, 2024)
historical_extremes <- history_long %>%
  filter(year %in% hist_years) %>%
  group_by(year, sex) %>%
  summarise(
    n_states = n_distinct(state),
    min_state = state[which.min(e0)],
    min_e0 = min(e0),
    max_state = state[which.max(e0)],
    max_e0 = max(e0),
    range_years = max(e0) - min(e0),
    mean_e0 = mean(e0),
    median_e0 = median(e0),
    sd_e0 = sd(e0),
    .groups = "drop"
  ) %>%
  arrange(year, sex)
write_table(historical_extremes, "historical_extremes_and_ranges_1971_2024.csv")

historical_by_group <- history_long %>%
  filter(year %in% hist_years) %>%
  group_by(year, sex, health_index_group) %>%
  summarise(
    n_states = n_distinct(state),
    mean_e0 = mean(e0),
    sd_e0 = sd(e0),
    median_e0 = median(e0),
    min_state = state[which.min(e0)],
    min_e0 = min(e0),
    max_state = state[which.max(e0)],
    max_e0 = max(e0),
    range_years = max(e0) - min(e0),
    .groups = "drop"
  ) %>%
  arrange(year, sex, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")))
write_table(historical_by_group, "historical_health_index_group_summary_1971_2024.csv")

# Use the wide state-year data here so both sexes are retained at each endpoint.
first_by_state <- df %>%
  group_by(state) %>%
  slice_min(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(state, first_year = year, first_Female = e0_f, first_Male = e0_m)
latest_by_state_wide <- df %>%
  group_by(state) %>%
  slice_max(year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(state, latest_year = year, latest_Female = e0_f, latest_Male = e0_m)

state_sex_gaps <- full_join(first_by_state, latest_by_state_wide, by = "state") %>%
  left_join(state_group_lookup %>% select(state, health_index_group), by = "state") %>%
  mutate(
    first_female_minus_male = first_Female - first_Male,
    latest_female_minus_male = latest_Female - latest_Male
  ) %>%
  select(state, health_index_group, first_year, first_Female, first_Male,
         first_female_minus_male, latest_year, latest_Female, latest_Male,
         latest_female_minus_male) %>%
  arrange(factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")), state)
write_table(state_sex_gaps, "historical_first_latest_values_and_sex_gaps_by_state.csv")

sex_gap_counts <- bind_rows(
  state_sex_gaps %>% transmute(endpoint = "Earliest observation", state, gap = first_female_minus_male),
  state_sex_gaps %>% transmute(endpoint = "Latest available observation", state, gap = latest_female_minus_male)
) %>%
  group_by(endpoint) %>%
  summarise(
    n_states = n(),
    male_life_expectancy_higher = sum(gap < 0, na.rm = TRUE),
    female_life_expectancy_higher = sum(gap > 0, na.rm = TRUE),
    equal = sum(gap == 0, na.rm = TRUE),
    narrowest_positive_female_advantage = if (any(gap > 0, na.rm = TRUE)) min(gap[gap > 0], na.rm = TRUE) else NA_real_,
    widest_positive_female_advantage = if (any(gap > 0, na.rm = TRUE)) max(gap[gap > 0], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )
write_table(sex_gap_counts, "historical_sex_gap_counts.csv")

# -----------------------------------------------------------------------------
# 4. Results Section 2: projections, dispersion, sex gaps, thresholds, equilibrium
# -----------------------------------------------------------------------------

forecast <- read_required("tables/main_regularised/forecast_summary.csv") %>%
  mutate(
    state = as.character(state),
    sex = as.character(sex),
    year = as.numeric(year),
    forecast_type = as.character(forecast_type)
  )

latent_forecast <- forecast %>%
  filter(forecast_type == "Latent trajectory") %>%
  left_join(state_group_lookup %>% select(state, health_index_group), by = "state")

selected_horizon_forecasts <- latent_forecast %>%
  filter(year %in% c(2030, 2040, 2050)) %>%
  select(any_of(c("state", "health_index_group", "sex", "year", "median", "q01", "q025", "q05", "q25", "q75", "q95", "q975", "q99", "width50", "width90", "width95"))) %>%
  arrange(sex, year, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")), state)
write_table(selected_horizon_forecasts, "latent_projection_estimates_by_state_2030_2040_2050.csv")

projection_2050 <- latent_forecast %>%
  filter(year == 2050, is.finite(median)) %>%
  select(state, health_index_group, sex, year, median, q025, q975, any_of(c("q01", "q05", "q25", "q75", "q95", "q99", "width50", "width90", "width95"))) %>%
  arrange(sex, median)
if (nrow(projection_2050) == 0L) stop("No 2050 latent forecasts in forecast_summary.csv.", call. = FALSE)
write_table(projection_2050, "latent_projection_estimates_2050_by_state.csv")

projection_group_summary <- projection_2050 %>%
  group_by(health_index_group, sex) %>%
  summarise(
    n_states = n_distinct(state),
    mean_state_posterior_median = mean(median),
    median_of_state_posterior_medians = median(median),
    lowest_state = state[which.min(median)],
    lowest_state_posterior_median = min(median),
    highest_state = state[which.max(median)],
    highest_state_posterior_median = max(median),
    mean_state_95CrI_width = mean(q975 - q025),
    median_state_95CrI_width = median(q975 - q025),
    .groups = "drop"
  ) %>%
  arrange(sex, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")))
write_table(projection_group_summary, "projection_2050_summary_by_health_index_group.csv")

# Compare observed range among states with an actual 2024 observation against the
# range of state posterior medians in 2050. These state counts can differ.
range_2024 <- history_long %>%
  filter(year == 2024) %>%
  group_by(sex) %>%
  summarise(
    comparison_period = "2024 observed",
    n_states = n_distinct(state),
    min_state = state[which.min(e0)], min_e0 = min(e0),
    max_state = state[which.max(e0)], max_e0 = max(e0),
    interstate_range = max(e0) - min(e0),
    .groups = "drop"
  )
range_2050 <- projection_2050 %>%
  group_by(sex) %>%
  summarise(
    comparison_period = "2050 projected posterior medians",
    n_states = n_distinct(state),
    min_state = state[which.min(median)], min_e0 = min(median),
    max_state = state[which.max(median)], max_e0 = max(median),
    interstate_range = max(median) - min(median),
    .groups = "drop"
  )
interstate_dispersion <- bind_rows(range_2024, range_2050) %>%
  group_by(sex) %>%
  mutate(
    change_in_range_vs_2024 = interstate_range - interstate_range[comparison_period == "2024 observed"][1],
    percent_change_in_range_vs_2024 = 100 * change_in_range_vs_2024 / interstate_range[comparison_period == "2024 observed"][1]
  ) %>%
  ungroup()
write_table(interstate_dispersion, "interstate_dispersion_2024_observed_vs_2050_projected.csv")

female_2050 <- projection_2050 %>% filter(sex == "Female") %>% select(state, health_index_group, female_median = median)
male_2050 <- projection_2050 %>% filter(sex == "Male") %>% select(state, male_median = median)
sex_gap_2050 <- inner_join(female_2050, male_2050, by = "state") %>%
  mutate(female_minus_male = female_median - male_median) %>%
  arrange(female_minus_male)
write_table(sex_gap_2050, "projected_female_minus_male_gap_2050_by_state.csv")

sex_gap_2050_group <- sex_gap_2050 %>%
  group_by(health_index_group) %>%
  summarise(
    n_states = n(),
    mean_female_minus_male = mean(female_minus_male),
    median_female_minus_male = median(female_minus_male),
    smallest_gap_state = state[which.min(female_minus_male)],
    smallest_gap = min(female_minus_male),
    largest_gap_state = state[which.max(female_minus_male)],
    largest_gap = max(female_minus_male),
    n_states_gap_over_5_years = sum(female_minus_male > 5),
    .groups = "drop"
  )
write_table(sex_gap_2050_group, "projected_female_minus_male_gap_2050_by_health_index_group.csv")

# State-specific equilibrium target estimates and difference from the latest SRS value.
equilibrium <- read_required("tables/main_regularised/state_equilibrium_targets.csv") %>%
  mutate(state = as.character(state), sex = as.character(sex), tier = as.character(tier)) %>%
  left_join(state_group_lookup %>% select(state, health_index_group), by = "state")
latest_long <- latest_by_state_wide %>%
  pivot_longer(
    cols = c(latest_Female, latest_Male),
    names_to = "sex",
    values_to = "latest_observed_e0"
  ) %>%
  mutate(
    sex = recode(sex, latest_Female = "Female", latest_Male = "Male"),
    latest_observed_year = latest_year
  ) %>%
  select(state, latest_observed_year, sex, latest_observed_e0)

# Because SRS input uncertainty is not supplied to the likelihood, these gap
# intervals translate the target interval by a fixed latest-observed value.
equilibrium_targets_by_state <- equilibrium %>%
  inner_join(latest_long, by = c("state", "sex")) %>%
  mutate(
    target_minus_latest_observation = median - latest_observed_e0,
    gap_q025 = q025 - latest_observed_e0,
    gap_q975 = q975 - latest_observed_e0
  ) %>%
  select(state, health_index_group, tier, sex, median, q025, q975,
         latest_observed_year, latest_observed_e0,
         target_minus_latest_observation, gap_q025, gap_q975) %>%
  arrange(sex, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")), state)
write_table(equilibrium_targets_by_state, "state_equilibrium_targets_and_gaps.csv")

equilibrium_group_summary <- equilibrium_targets_by_state %>%
  group_by(health_index_group, sex) %>%
  summarise(
    n_states = n_distinct(state),
    mean_state_target_median = mean(median),
    median_state_target_median = median(median),
    lowest_target_state = state[which.min(median)],
    lowest_target_median = min(median),
    highest_target_state = state[which.max(median)],
    highest_target_median = max(median),
    mean_target_minus_latest_observation = mean(target_minus_latest_observation),
    median_target_minus_latest_observation = median(target_minus_latest_observation),
    .groups = "drop"
  )
write_table(equilibrium_group_summary, "equilibrium_target_summary_by_health_index_group.csv")

# Direct posterior draws of the three Health Index tier centres (transformed parameters).
main_extract_tier <- tryCatch(
  rstan::extract(main_fit, pars = c("C_tier_f", "C_tier_m", "C_f", "C_m", "delta_C", "kappa_f", "kappa_m", "lambda_dev", "sigma_obs", "sigma_dev", "tau_C"), permuted = TRUE),
  error = function(e) {
    message("Could not extract some optional tier/parameter draws: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(main_extract_tier) && !is.null(main_extract_tier$C_tier_f) && !is.null(main_extract_tier$C_tier_m)) {
  tier_names <- c("Front-runner", "Achiever", "Aspirant")
  tier_summary <- map_dfr(seq_along(tier_names), function(j) {
    f <- as.numeric(main_extract_tier$C_tier_f[, j])
    m <- as.numeric(main_extract_tier$C_tier_m[, j])
    bind_rows(
      tibble(health_index_group = tier_names[j], sex = "Female",
             median = median(f), q025 = unname(quantile(f, .025)), q975 = unname(quantile(f, .975))),
      tibble(health_index_group = tier_names[j], sex = "Male",
             median = median(m), q025 = unname(quantile(m, .025)), q975 = unname(quantile(m, .975)))
    )
  })
  write_table(tier_summary, "health_index_tier_equilibrium_centres.csv")
}

# Threshold probabilities: computed by Script 04 using the latent trajectory.
threshold_state <- read_required("tables/threshold_probabilities_2024_2050/state_probability_below_80_2024_2050.csv") %>%
  mutate(state = as.character(state), sex = as.character(sex), year = as.integer(year))
write_table(threshold_state, "threshold_probability_below_80_by_state_2024_2050.csv")

threshold_change <- threshold_state %>%
  select(state, health_index_group, sex, year, p_e0_below_80, p_e0_at_least_80, posterior_median, q025, q975) %>%
  pivot_wider(
    names_from = year,
    values_from = c(p_e0_below_80, p_e0_at_least_80, posterior_median, q025, q975),
    names_glue = "{.value}_{year}"
  ) %>%
  mutate(
    change_p_below_80 = p_e0_below_80_2050 - p_e0_below_80_2024,
    change_p_at_least_80 = p_e0_at_least_80_2050 - p_e0_at_least_80_2024
  ) %>%
  arrange(sex, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")), state)
write_table(threshold_change, "threshold_probability_change_2024_to_2050_by_state.csv")

threshold_group <- read_required("tables/threshold_probabilities_2024_2050/health_index_group_probability_below_80_2024_2050.csv")
write_table(threshold_group, "threshold_probability_below_80_by_health_index_group.csv")
threshold_group_change <- threshold_group %>%
  mutate(year = as.integer(year)) %>%
  select(health_index_group, sex, year, mean_state_probability_below_80,
         median_state_probability_below_80, n_states_probability_below_80_over_50,
         pct_states_probability_below_80_over_50) %>%
  pivot_wider(
    names_from = year,
    values_from = c(mean_state_probability_below_80, median_state_probability_below_80,
                    n_states_probability_below_80_over_50, pct_states_probability_below_80_over_50),
    names_glue = "{.value}_{year}"
  ) %>%
  mutate(
    change_mean_state_probability = mean_state_probability_below_80_2050 - mean_state_probability_below_80_2024,
    change_median_state_probability = median_state_probability_below_80_2050 - median_state_probability_below_80_2024
  )
write_table(threshold_group_change, "threshold_probability_change_by_health_index_group.csv")

# -----------------------------------------------------------------------------
# 5. Results Section 3: COVID shock and post-2021 decay
# -----------------------------------------------------------------------------

covid_vars <- c("delta2020_f", "delta2021_f", "delta2020_m", "delta2021_m", "phi_covid_f", "phi_covid_m")
covid_extract <- rstan::extract(main_fit, pars = covid_vars, permuted = TRUE)

summarise_draw_vector <- function(x) {
  x <- as.numeric(x)
  tibble(
    median = median(x),
    q025 = unname(quantile(x, .025, names = FALSE, type = 8)),
    q975 = unname(quantile(x, .975, names = FALSE, type = 8)),
    mean = mean(x),
    sd = sd(x)
  )
}

covid_rows <- list()
add_covid_summary <- function(draws, quantity, sex, unit) {
  s <- summarise_draw_vector(draws)
  s %>% mutate(quantity = quantity, sex = sex, unit = unit, .before = 1)
}
covid_shock_summary <- bind_rows(
  add_covid_summary(covid_extract$delta2020_f, "2020 deviation", "Female", "years"),
  add_covid_summary(covid_extract$delta2021_f, "2021 deviation", "Female", "years"),
  add_covid_summary(covid_extract$delta2020_m, "2020 deviation", "Male", "years"),
  add_covid_summary(covid_extract$delta2021_m, "2021 deviation", "Male", "years"),
  add_covid_summary(covid_extract$phi_covid_f, "Post-2021 annual decay fraction", "Female", "proportion"),
  add_covid_summary(covid_extract$phi_covid_m, "Post-2021 annual decay fraction", "Male", "proportion"),
  add_covid_summary(log(0.5) / log(covid_extract$phi_covid_f), "Post-2021 shock half-life", "Female", "years"),
  add_covid_summary(log(0.5) / log(covid_extract$phi_covid_m), "Post-2021 shock half-life", "Male", "years")
)
write_table(covid_shock_summary, "covid_shock_and_decay_posterior_summary.csv")

# Export selected model parameters as raw posterior summaries from the saved
# RStan summary, including R-hat and ESS columns where available.
selected_parameter_pattern <- paste0(
  "^(C_f|C_m|delta_C|C_tier_f|C_tier_m|C_state_f|C_state_m|tau_C|",
  "kappa_f|kappa_m|lambda_dev|deviation_half_life|convergence_half_life_f|",
  "convergence_half_life_m|years_to_90pct_f|years_to_90pct_m|sigma_obs|sigma_dev|",
  "long_run_sd_dev_f|long_run_sd_dev_m|delta2020_f|delta2021_f|delta2020_m|",
  "delta2021_m|phi_covid_f|phi_covid_m)(\\[|$)"
)
main_mcmc <- read_required("tables/main_regularised/mcmc_summary_all_quantities.csv", name_repair = "minimal")
selected_params <- main_mcmc %>%
  filter(grepl(selected_parameter_pattern, parameter)) %>%
  arrange(parameter)
write_table(selected_params, "selected_main_model_parameter_summaries.csv")

# -----------------------------------------------------------------------------
# 6. Results Section 4: out-of-sample predictive performance
# -----------------------------------------------------------------------------

holdout_overall <- read_required("tables/holdout_validation/holdout_metrics_overall.csv") %>%
  mutate(forecast_start_year = as.integer(forecast_start_year)) %>%
  arrange(forecast_start_year, sex, model)
write_table(holdout_overall, "holdout_metrics_overall_all_models.csv")

holdout_by_lead <- read_required("tables/holdout_validation/holdout_metrics_by_origin_sex_lead.csv") %>%
  mutate(forecast_start_year = as.integer(forecast_start_year)) %>%
  arrange(forecast_start_year, sex, model, lead_band)
write_table(holdout_by_lead, "holdout_metrics_by_origin_sex_lead.csv")

holdout_by_state <- read_optional("tables/holdout_validation/holdout_metrics_by_state.csv")
if (!is.null(holdout_by_state)) write_table(holdout_by_state, "holdout_metrics_by_state.csv")
holdout_predictions <- read_optional("tables/holdout_validation/holdout_predictions_main_model_and_benchmarks.csv")
if (!is.null(holdout_predictions)) write_table(holdout_predictions, "holdout_predictions_all_models_raw.csv")
holdout_diag <- read_optional("tables/holdout_validation/sampler_diagnostics_main_model_holdout_fits.csv")
if (!is.null(holdout_diag)) write_table(holdout_diag, "holdout_sampler_diagnostics_by_origin.csv")
excluded_states <- read_optional("tables/holdout_validation/states_without_pre_origin_training_data.csv")
if (!is.null(excluded_states)) write_table(excluded_states, "holdout_states_without_pre_origin_training_data.csv")

# -----------------------------------------------------------------------------
# 7. Results Section 5: diagnostics and predictive checks
# -----------------------------------------------------------------------------

# Combine full-data posterior sampler and generated-quantity checks.
posterior_diag_rows <- list()
gq_diag_rows <- list()
posterior_parameter_diag_rows <- list()
for (d in names(model_dirs)) {
  label <- unname(model_dirs[d])
  diag_rel <- file.path("tables", d, "diagnostic_status.csv")
  gq_rel <- file.path("tables", d, "generated_quantities_status.csv")
  mcmc_rel <- file.path("tables", d, "mcmc_summary_all_quantities.csv")

  diag_dat <- read_optional(diag_rel)
  if (!is.null(diag_dat)) posterior_diag_rows[[length(posterior_diag_rows) + 1L]] <- diag_dat %>% mutate(model = label, .before = 1)

  gq_dat <- read_optional(gq_rel)
  if (!is.null(gq_dat)) gq_diag_rows[[length(gq_diag_rows) + 1L]] <- gq_dat %>% mutate(model = label, .before = 1)

  mcmc_dat <- read_optional(mcmc_rel, name_repair = "minimal")
  if (!is.null(mcmc_dat) && all(c("parameter", "Rhat", "n_eff") %in% names(mcmc_dat))) {
    # Same exclusion pattern used by the workflow's diagnostic screen. This keeps
    # very large generated quantities out of the R-hat/ESS parameter diagnostics.
    omit_pattern <- paste0(
      "^(y_rep_f|y_rep_m|mu_curve_f|mu_curve_m|C_tier_f|C_tier_m|tier_effect_f|tier_effect_m|",
      "elapsed_since_first|baseline_f|baseline_m|covid_shock_f|covid_shock_m|",
      "trend_forecast_f|trend_forecast_m|mean_forecast_f|mean_forecast_m|",
      "latent_forecast_f|latent_forecast_m|e0_forecast_f|e0_forecast_m|",
      "u_pred_mean|u_pred_cov|u_filter_mean|u_filter_cov|marginal_loglik|",
      "I2|L_scaled|L_Omega_dev|L_Omega_C|Omega_dev|Omega_C|L_C_scaled)(\\[|$)"
    )
    diag_params <- mcmc_dat %>%
      filter(parameter != "lp__") %>%
      filter(!grepl(omit_pattern, parameter) | parameter %in% c("L_Omega_dev[2,1]", "L_Omega_C[2,1]")) %>%
      mutate(model = label, .before = 1)
    posterior_parameter_diag_rows[[length(posterior_parameter_diag_rows) + 1L]] <- diag_params
  }
}
if (length(posterior_diag_rows)) write_table(bind_rows(posterior_diag_rows), "posterior_sampler_diagnostic_status_by_model.csv")
if (length(gq_diag_rows)) write_table(bind_rows(gq_diag_rows), "generated_quantities_status_by_model.csv")
if (length(posterior_parameter_diag_rows)) {
  par_diag_all <- bind_rows(posterior_parameter_diag_rows)
  write_table(par_diag_all, "posterior_parameter_rhat_ess_all_models.csv")
  par_diag_summary <- par_diag_all %>%
    group_by(model) %>%
    summarise(
      n_parameter_rows = n(),
      maximum_rhat = if (any(is.finite(Rhat))) max(Rhat[is.finite(Rhat)]) else NA_real_,
      minimum_ess = if (any(is.finite(n_eff))) min(n_eff[is.finite(n_eff)]) else NA_real_,
      n_rhat_above_1_01 = sum(Rhat > 1.01, na.rm = TRUE),
      n_ess_below_100 = sum(n_eff < 100, na.rm = TRUE),
      n_nonfinite_rhat = sum(!is.finite(Rhat)),
      n_nonfinite_ess = sum(!is.finite(n_eff)),
      .groups = "drop"
    )
  write_table(par_diag_summary, "posterior_parameter_rhat_ess_summary_by_model.csv")
}

# Prior-only diagnostic files from the three prior predictive checks.
prior_diag_rows <- list()
prior_gq_rows <- list()
for (pm in 2:4) {
  prior_label <- c(`2` = "Main prior", `3` = "Tighter prior", `4` = "Wide prior")[as.character(pm)]
  rel_diag <- file.path("tables", "prior_predictive", paste0("mode_", pm), "diagnostic_status.csv")
  rel_gq <- file.path("tables", "prior_predictive", paste0("mode_", pm), "generated_quantities_status.csv")
  pd <- read_optional(rel_diag)
  if (!is.null(pd)) prior_diag_rows[[length(prior_diag_rows) + 1L]] <- pd %>% mutate(prior_setting = unname(prior_label), prior_mode = pm, .before = 1)
  pg <- read_optional(rel_gq)
  if (!is.null(pg)) prior_gq_rows[[length(prior_gq_rows) + 1L]] <- pg %>% mutate(prior_setting = unname(prior_label), prior_mode = pm, .before = 1)
}
if (length(prior_diag_rows)) write_table(bind_rows(prior_diag_rows), "prior_predictive_sampler_diagnostic_status.csv")
if (length(prior_gq_rows)) write_table(bind_rows(prior_gq_rows), "prior_predictive_generated_quantities_status.csv")

# Overall posterior predictive coverage table already emitted by the engine.
ppc_overall <- read_required("tables/main_regularised/historical_one_step_ppc_coverage.csv")
write_table(ppc_overall, "posterior_predictive_coverage_overall_by_sex.csv")

# Compute state-by-sex PPC coverage from stored y_rep draws. Initial observations
# are excluded because they anchor the trend and y_rep is set to the observed value.
ppc_by_state <- NULL
if ("is_first" %in% names(stan_data)) {
  ppc_idx <- which(as.integer(stan_data$is_first) == 0L)
  ppc_ex <- rstan::extract(main_fit, pars = c("y_rep_f", "y_rep_m"), permuted = TRUE)
  ppc_rows <- list()
  for (sx in c("Female", "Male")) {
    suffix <- if (sx == "Female") "f" else "m"
    reps <- ppc_ex[[paste0("y_rep_", suffix)]]
    if (length(dim(reps)) != 2L || dim(reps)[2] != nrow(df)) {
      warning("Could not calculate state-specific PPC coverage: unexpected y_rep dimensions.")
      ppc_rows <- list()
      break
    }
    obs <- if (sx == "Female") df$e0_f else df$e0_m
    q <- t(apply(reps[, ppc_idx, drop = FALSE], 2, quantile,
                 probs = c(.05, .25, .50, .75, .95), names = FALSE, type = 8))
    ppc_rows[[length(ppc_rows) + 1L]] <- tibble(
      state = df$state[ppc_idx],
      sex = sx,
      observed = obs[ppc_idx],
      q05 = q[, 1], q25 = q[, 2], q50 = q[, 3], q75 = q[, 4], q95 = q[, 5],
      covered50 = obs[ppc_idx] >= q[, 2] & obs[ppc_idx] <= q[, 4],
      covered90 = obs[ppc_idx] >= q[, 1] & obs[ppc_idx] <= q[, 5]
    )
  }
  if (length(ppc_rows)) {
    ppc_by_state <- bind_rows(ppc_rows) %>%
      group_by(state, sex) %>%
      summarise(
        n_observations = n(),
        coverage50 = mean(covered50),
        coverage90 = mean(covered90),
        mean_observed_e0 = mean(observed),
        mean_predictive_median = mean(q50),
        .groups = "drop"
      ) %>%
      left_join(state_group_lookup %>% select(state, health_index_group), by = "state") %>%
      arrange(sex, factor(health_index_group, levels = c("Front-runner", "Achiever", "Aspirant")), state)
    write_table(ppc_by_state, "posterior_predictive_coverage_by_state_and_sex.csv")
  }
}

prior_ppc <- read_optional("tables/prior_predictive/prior_predictive_historical_interval_coverage.csv")
if (!is.null(prior_ppc)) write_table(prior_ppc, "prior_predictive_historical_coverage_by_prior_and_sex.csv")
prior_plausibility <- read_optional("tables/prior_predictive/prior_predictive_plausibility_shares.csv")
if (!is.null(prior_plausibility)) {
  prior_plausibility_selected <- prior_plausibility %>% filter(year %in% c(2030, 2040, 2050))
  write_table(prior_plausibility_selected, "prior_predictive_plausibility_shares_2030_2040_2050.csv")
}
prior_ranges <- read_optional("tables/prior_predictive/prior_predictive_quantile_ranges_001_099.csv")
if (!is.null(prior_ranges)) write_table(prior_ranges, "prior_predictive_quantile_ranges.csv")

# -----------------------------------------------------------------------------
# 8. Prior sensitivity results
# -----------------------------------------------------------------------------

sensitivity_compare <- read_optional("tables/prior_sensitivity_comparison/state_year_numeric_sensitivity_comparison.csv")
if (!is.null(sensitivity_compare)) {
  write_table(sensitivity_compare, "prior_sensitivity_state_year_comparison_all_horizons.csv")
  sensitivity_2050 <- sensitivity_compare %>% filter(year == 2050) %>%
    select(any_of(c("forecast_type", "sex", "state", "g", "year",
                    "median_mode2", "q025_mode2", "q975_mode2",
                    "median_mode3", "q025_mode3", "q975_mode3",
                    "median_mode4", "q025_mode4", "q975_mode4",
                    "shift_tight_vs_main", "shift_wide_vs_main",
                    "width95_change_tight", "width95_change_wide",
                    "interval95_overlap_tight", "interval95_overlap_wide")))
  write_table(sensitivity_2050, "prior_sensitivity_comparison_2050.csv")

  sensitivity_shift_2050 <- sensitivity_2050 %>%
    group_by(forecast_type, sex) %>%
    summarise(
      n_states = n_distinct(state),
      median_absolute_shift_tighter_prior = median(abs(shift_tight_vs_main), na.rm = TRUE),
      mean_absolute_shift_tighter_prior = mean(abs(shift_tight_vs_main), na.rm = TRUE),
      maximum_absolute_shift_tighter_prior = max(abs(shift_tight_vs_main), na.rm = TRUE),
      n_states_tighter_shift_over_1_year = sum(abs(shift_tight_vs_main) > 1, na.rm = TRUE),
      median_absolute_shift_wider_prior = median(abs(shift_wide_vs_main), na.rm = TRUE),
      mean_absolute_shift_wider_prior = mean(abs(shift_wide_vs_main), na.rm = TRUE),
      maximum_absolute_shift_wider_prior = max(abs(shift_wide_vs_main), na.rm = TRUE),
      n_states_wider_shift_over_1_year = sum(abs(shift_wide_vs_main) > 1, na.rm = TRUE),
      n_state_intervals_overlap_tighter = sum(interval95_overlap_tight, na.rm = TRUE),
      n_state_intervals_overlap_wider = sum(interval95_overlap_wide, na.rm = TRUE),
      .groups = "drop"
    )
  write_table(sensitivity_shift_2050, "prior_sensitivity_2050_shift_summary.csv")
}

sens_width <- read_optional("tables/prior_sensitivity_comparison/summary_interval_widths_by_horizon.csv")
if (!is.null(sens_width)) write_table(sens_width, "prior_sensitivity_interval_widths_by_horizon.csv")
sens_shift <- read_optional("tables/prior_sensitivity_comparison/summary_median_shifts_by_horizon.csv")
if (!is.null(sens_shift)) write_table(sens_shift, "prior_sensitivity_median_shifts_by_horizon.csv")
sens_ppc <- read_optional("tables/prior_sensitivity_comparison/posterior_predictive_coverage_by_prior_and_sex.csv")
if (!is.null(sens_ppc)) write_table(sens_ppc, "prior_sensitivity_posterior_predictive_coverage.csv")

# -----------------------------------------------------------------------------
# 9. A compact table of headline numeric results for quick review
# -----------------------------------------------------------------------------

key_rows <- list()
add_key <- function(section, quantity, value, sex = NA_character_, year = NA_real_, group = NA_character_, unit = "", note = "") {
  key_rows[[length(key_rows) + 1L]] <<- tibble(
    section = section, quantity = quantity, sex = sex, year = year,
    health_index_group = group, value = as.numeric(value), unit = unit, note = note
  )
}

for (i in seq_len(nrow(historical_extremes))) {
  r <- historical_extremes[i, ]
  add_key("Historical", "Number of states", sex = r$sex, year = r$year, value = r$n_states, unit = "states")
  add_key("Historical", "Minimum e0", sex = r$sex, year = r$year, value = r$min_e0, unit = "years", note = r$min_state)
  add_key("Historical", "Maximum e0", sex = r$sex, year = r$year, value = r$max_e0, unit = "years", note = r$max_state)
  add_key("Historical", "Interstate range", sex = r$sex, year = r$year, value = r$range_years, unit = "years")
}
for (i in seq_len(nrow(interstate_dispersion))) {
  r <- interstate_dispersion[i, ]
  add_key("Projections", paste0("Interstate range: ", r$comparison_period), sex = r$sex, year = ifelse(r$comparison_period == "2024 observed", 2024, 2050), value = r$interstate_range, unit = "years", note = paste0(r$min_state, " to ", r$max_state, "; n=", r$n_states))
}
for (i in seq_len(nrow(projection_group_summary))) {
  r <- projection_group_summary[i, ]
  add_key("Projections", "Mean of state posterior medians in 2050", value = r$mean_state_posterior_median,
          sex = r$sex, year = 2050, group = r$health_index_group, unit = "years",
          note = paste0("State medians range from ", r$lowest_state, " (", signif(r$lowest_state_posterior_median, 4), ") to ", r$highest_state, " (", signif(r$highest_state_posterior_median, 4), ")"))
}
for (i in seq_len(nrow(threshold_group))) {
  r <- threshold_group[i, ]
  add_key("Projections / threshold probabilities", "Mean state-specific P(e0 < 80)", value = r$mean_state_probability_below_80,
          sex = r$sex, year = as.integer(r$year), group = r$health_index_group, unit = "probability",
          note = "Equal-weighted average across states; not a population-weighted group probability")
  add_key("Projections / threshold probabilities", "States with P(e0 < 80) > 0.5", value = r$n_states_probability_below_80_over_50,
          sex = r$sex, year = as.integer(r$year), group = r$health_index_group, unit = "states")
}
for (i in seq_len(nrow(sex_gap_2050_group))) {
  r <- sex_gap_2050_group[i, ]
  add_key("Projections / sex gap", "Median state female-minus-male gap", value = r$median_female_minus_male,
          year = 2050, group = r$health_index_group, unit = "years")
}
for (i in seq_len(nrow(covid_shock_summary))) {
  r <- covid_shock_summary[i, ]
  add_key("COVID shock", r$quantity, sex = r$sex, value = r$median, unit = r$unit,
          note = paste0("95% CrI ", signif(r$q025, 4), " to ", signif(r$q975, 4)))
}
for (i in seq_len(nrow(holdout_overall))) {
  r <- holdout_overall[i, ]
  for (metric in intersect(c("MAE", "MSE", "RMSE", "CRPS", "coverage90"), names(holdout_overall))) {
    add_key("Out-of-sample validation", metric, sex = r$sex, year = r$forecast_start_year,
            value = r[[metric]], unit = ifelse(metric == "MSE", "years squared", ifelse(metric == "coverage90", "proportion", "years")),
            note = as.character(r$model))
  }
}
write_table(bind_rows(key_rows), "key_numbers_for_results_review.csv")

# -----------------------------------------------------------------------------
# 10. README and status summary
# -----------------------------------------------------------------------------

status_lines <- c(
  "MANUSCRIPT RESULTS RAW-NUMBER EXTRACTION",
  paste0("Project root: ", project_root),
  paste0("Output directory: ", out_dir),
  "No posterior model was fitted or refitted by this script.",
  "",
  "Main output files:",
  "- Historical: historical_extremes_and_ranges_1971_2024.csv; historical_health_index_group_summary_1971_2024.csv; historical_state_values_1971_2024.csv; historical_first_latest_values_and_sex_gaps_by_state.csv; historical_sex_gap_counts.csv",
  "- Projections: latent_projection_estimates_by_state_2030_2040_2050.csv; latent_projection_estimates_2050_by_state.csv; projection_2050_summary_by_health_index_group.csv; interstate_dispersion_2024_observed_vs_2050_projected.csv; projected_female_minus_male_gap_2050_by_state.csv",
  "- Equilibrium: state_equilibrium_targets_and_gaps.csv; equilibrium_target_summary_by_health_index_group.csv; health_index_tier_equilibrium_centres.csv when available",
  "- Threshold probability: threshold_probability_below_80_by_state_2024_2050.csv; threshold_probability_change_2024_to_2050_by_state.csv; threshold_probability_below_80_by_health_index_group.csv; threshold_probability_change_by_health_index_group.csv",
  "- COVID: covid_shock_and_decay_posterior_summary.csv",
  "- Validation: holdout_metrics_overall_all_models.csv; holdout_metrics_by_origin_sex_lead.csv; holdout_metrics_by_state.csv when available; holdout predictions and diagnostics when available",
  "- Diagnostics/PPC: posterior_sampler_diagnostic_status_by_model.csv; generated_quantities_status_by_model.csv; posterior_parameter_rhat_ess_all_models.csv; posterior_predictive_coverage_overall_by_sex.csv; posterior_predictive_coverage_by_state_and_sex.csv when available",
  "- Prior checks: prior predictive files and diagnostic status, where available",
  "- Sensitivity: prior_sensitivity_comparison_2050.csv and summary tables, where available",
  "- key_numbers_for_results_review.csv is a compact roll-up. Consult the fuller tables for complete raw estimates.",
  "",
  "Definitions and cautions:",
  "1. Forecast projections and threshold probabilities refer to the latent trajectory, not future observation draws.",
  "2. 2024 threshold probabilities for states lacking a 2024 observation are model-based latent-state estimates propagated from the latest observation.",
  "3. Health Index group means are state-level/equal-weighted summaries, not population-weighted pooled life expectancy.",
  "4. 2024 observed interstate ranges use only states with actual 2024 observations; 2050 ranges use states with 2050 projections. Check n_states when comparing ranges.",
  "5. Target-minus-latest-observation intervals translate uncertainty in the target while treating the latest observed SRS value as fixed; input SRS standard errors are not modelled.",
  "6. Point-only holdout benchmarks have no CRPS or predictive interval coverage; these cells should remain NA.",
  "7. Convergence diagnostics should be interpreted separately from predictive validation and posterior predictive checks.",
  "8. Missing optional files are listed as exists=FALSE in input_file_inventory.csv; do not interpret missing outputs as passed checks."
)
writeLines(status_lines, con = file.path(out_dir, "README_results_raw.txt"))

cat("\nRESULTS EXTRACTION COMPLETE\n")
cat("Project root: ", project_root, "\n", sep = "")
cat("Raw tables:   ", out_dir, "\n", sep = "")
cat("No refitting was performed. Review input_file_inventory.csv for any missing optional inputs.\n")
