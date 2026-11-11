# 04_threshold_probabilities_2050.R (2024 and 2050 heatmap)
#
# State-specific posterior probability that life expectancy is below 80 years
# in 2024 and 2050, shown in a heatmap faceted by sex.
#
# Primary quantity shown:
#   P(e0 < 80)
# The complementary probability P(e0 >= 80) = 1 - P(e0 < 80) is also saved
# in the output tables, but is not plotted separately because it is redundant.
#
# Uses the latent trajectory (trend + latent OU deviation), not the future-
# observation distribution that adds observation-level residual noise.
#
# 2024 handling:
#   - If a state's latest observation is in 2024, sample the latent state from
#     its filtered posterior conditional on that state's observed history.
#   - If its latest observation predates 2024, propagate the latent OU state
#     from its latest observation to 2024 using the fitted model.
# This script assumes the input panel has no observations after 2024.
#
# Inputs:
#   fits/prep_objects.rds
#   fits/fit_main_regularised.rds
#
# Outputs:
#   figures/threshold_probabilities_2024_2050/probability_below_80_heatmap_2024_2050.png
#   tables/threshold_probabilities_2024_2050/state_probability_below_80_2024_2050.csv
#   tables/threshold_probabilities_2024_2050/health_index_group_probability_below_80_2024_2050.csv
#
# Run from the project root:
#   Rscript codes/04_threshold_probabilities_2050.R

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
})

set.seed(20261012)

# -----------------------------------------------------------------------------
# 1. Paths and input checks
# -----------------------------------------------------------------------------

project_root <- here::here()
prep_path <- file.path(project_root, "fits", "prep_objects.rds")
fit_path <- file.path(project_root, "fits", "fit_main_regularised.rds")
figure_dir <- file.path(project_root, "figures", "threshold_probabilities_2024_2050")
table_dir <- file.path(project_root, "tables", "threshold_probabilities_2024_2050")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(prep_path)) stop("Preparation object not found: ", prep_path, call. = FALSE)
if (!file.exists(fit_path)) stop("Main-model fit not found: ", fit_path, call. = FALSE)

prep <- readRDS(prep_path)
fit <- readRDS(fit_path)
stan_data <- prep$stan_data_full
df <- tibble::as_tibble(prep$df)
state_names <- as.character(prep$state_levels_full)
N_states <- length(state_names)

if (nrow(df) != stan_data$N) stop("prep$df rows do not match Stan data N.", call. = FALSE)
if (length(stan_data$last_idx) != N_states || length(stan_data$last_year) != N_states) {
  stop("State ordering/index mismatch between preparation object and Stan data.", call. = FALSE)
}

if (any(stan_data$last_year > 2024)) {
  bad <- state_names[stan_data$last_year > 2024]
  stop(
    "This script estimates 2024 using each state's final available observation, so it assumes no observations after 2024. Found later observations for: ",
    paste(bad, collapse = ", "), call. = FALSE
  )
}

forecast_years <- as.numeric(stan_data$forecast_years)
idx_2050 <- match(2050, forecast_years)
if (is.na(idx_2050)) stop("2050 is not present in stan_data$forecast_years.", call. = FALSE)

# -----------------------------------------------------------------------------
# 2. Health Index group map: kept consistent with the model engine
# -----------------------------------------------------------------------------

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
state_key <- normalise_state(state_names)
group_idx <- match(state_key, state_tier_map$state_key)
if (anyNA(group_idx)) {
  stop("No Health Index group mapped for: ", paste(state_names[is.na(group_idx)], collapse = ", "), call. = FALSE)
}
state_groups <- state_tier_map$health_index_group[group_idx]

# Make sure all group labels in the plotting map match the states in the data.
state_order <- tibble(state = state_names, health_index_group = state_groups) %>%
  mutate(health_index_group = factor(
    health_index_group,
    levels = c("Front-runner", "Achiever", "Aspirant")
  )) %>%
  arrange(health_index_group, state) %>%
  pull(state)

# -----------------------------------------------------------------------------
# 3. Extract posterior quantities needed for 2024 and 2050
# -----------------------------------------------------------------------------

pars_needed <- c(
  "latent_forecast_f", "latent_forecast_m",
  "C_state_f", "C_state_m", "baseline_f", "baseline_m",
  "kappa_f", "kappa_m", "lambda_dev", "sigma_dev", "Omega_dev",
  "delta2021_f", "delta2021_m", "phi_covid_f", "phi_covid_m",
  "u_filter_mean", "u_filter_cov"
)
ex <- rstan::extract(fit, pars = pars_needed, permuted = TRUE)

required_extracted <- setdiff(pars_needed, names(ex))
if (length(required_extracted)) {
  stop(
    "The saved fit does not contain the required extracted quantities: ",
    paste(required_extracted, collapse = ", "),
    ". Check that this fit was generated from the current Stan model.", call. = FALSE
  )
}

latent_f_2050 <- ex$latent_forecast_f[, , idx_2050, drop = FALSE]
latent_m_2050 <- ex$latent_forecast_m[, , idx_2050, drop = FALSE]
dim(latent_f_2050) <- c(dim(ex$latent_forecast_f)[1], N_states)
dim(latent_m_2050) <- c(dim(ex$latent_forecast_m)[1], N_states)

n_draws <- nrow(latent_f_2050)
if (nrow(latent_m_2050) != n_draws) stop("Female/male 2050 draw counts differ.", call. = FALSE)
cat("Project root: ", project_root, "\n", sep = "")
cat("States: ", N_states, " | posterior draws per sex/year: ", n_draws, "\n", sep = "")
cat("Latest observed year across states: ", max(stan_data$last_year), "\n", sep = "")

# The filtering distribution gives p(u_t | observations through t, parameters).
# Since 2024 is the end of the observed panel, this conditions on the complete
# available state-specific series wherever the final observed year is 2024.
# For states with an earlier final observation, we propagate that state forward.

rmvn2 <- function(mu, cov) {
  cov <- (cov + t(cov)) / 2
  eg <- eigen(cov, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  as.numeric(mu + eg$vectors %*% (sqrt(vals) * rnorm(2)))
}

ou_diffusion_multiplier <- function(rate, delta_t) {
  if (rate * delta_t < 1e-8) return(sqrt(delta_t))
  sqrt(max(0, -expm1(-2 * rate * delta_t) / (2 * rate)))
}

# Return one posterior draw of the two-sex latent state at 2024 for each state.
# Uses the filtered state distribution at the last observed row, then propagates
# the continuous-time OU process if the state's last observation predates 2024.
latent_2024_f <- matrix(NA_real_, nrow = n_draws, ncol = N_states)
latent_2024_m <- matrix(NA_real_, nrow = n_draws, ncol = N_states)

for (d in seq_len(n_draws)) {
  for (gi in seq_len(N_states)) {
    last_i <- as.integer(stan_data$last_idx[gi])
    last_year <- as.numeric(stan_data$last_year[gi])
    first_year <- as.numeric(stan_data$first_year[gi])
    gap <- 2024 - last_year
    
    # Posterior draw of latent OU state at the last observed time.
    u_mean <- as.numeric(ex$u_filter_mean[d, last_i, ])
    u_cov <- matrix(ex$u_filter_cov[d, last_i, , ], nrow = 2, ncol = 2)
    u_last <- rmvn2(u_mean, u_cov)
    
    if (gap > 0) {
      lambda <- as.numeric(ex$lambda_dev[d])
      sigma <- as.numeric(ex$sigma_dev[d, ])
      omega <- matrix(ex$Omega_dev[d, , ], nrow = 2, ncol = 2)
      D <- diag(sigma, nrow = 2, ncol = 2)
      dt_mult <- ou_diffusion_multiplier(lambda, gap)
      Q <- (dt_mult^2) * D %*% omega %*% D
      u_2024 <- exp(-lambda * gap) * u_last + rmvn2(c(0, 0), Q)
    } else {
      u_2024 <- u_last
    }
    
    # Smooth trend plus the model's post-2021 COVID-shock decay at 2024.
    elapsed <- max(2024 - first_year, 0)
    covid_f <- as.numeric(ex$delta2021_f[d]) * as.numeric(ex$phi_covid_f[d])^3
    covid_m <- as.numeric(ex$delta2021_m[d]) * as.numeric(ex$phi_covid_m[d])^3
    
    trend_f <- as.numeric(ex$C_state_f[d, gi]) +
      (as.numeric(ex$baseline_f[d, gi]) - as.numeric(ex$C_state_f[d, gi])) *
      exp(-as.numeric(ex$kappa_f[d]) * elapsed) + covid_f
    
    trend_m <- as.numeric(ex$C_state_m[d, gi]) +
      (as.numeric(ex$baseline_m[d, gi]) - as.numeric(ex$C_state_m[d, gi])) *
      exp(-as.numeric(ex$kappa_m[d]) * elapsed) + covid_m
    
    latent_2024_f[d, gi] <- trend_f + u_2024[1]
    latent_2024_m[d, gi] <- trend_m + u_2024[2]
  }
}

if (!all(is.finite(latent_2024_f)) || !all(is.finite(latent_2024_m))) {
  stop("Non-finite 2024 latent draws were generated; refusing to produce the heatmap.", call. = FALSE)
}

# 2050 latent draws are taken directly from generated quantities of the main fit.
latent_2050_f <- matrix(latent_f_2050, nrow = n_draws, ncol = N_states)
latent_2050_m <- matrix(latent_m_2050, nrow = n_draws, ncol = N_states)

# -----------------------------------------------------------------------------
# 4. Compute probability of life expectancy below 80 years by state, sex, and year
# -----------------------------------------------------------------------------

probability_rows <- list()
k <- 0L
for (yr in c(2024L, 2050L)) {
  for (sx in c("Female", "Male")) {
    dm <- if (yr == 2024L && sx == "Female") latent_2024_f else
      if (yr == 2024L && sx == "Male") latent_2024_m else
        if (yr == 2050L && sx == "Female") latent_2050_f else latent_2050_m
    
    for (gi in seq_len(N_states)) {
      z <- dm[, gi]
      k <- k + 1L
      probability_rows[[k]] <- tibble(
        state = state_names[gi],
        health_index_group = state_groups[gi],
        sex = sx,
        year = yr,
        p_e0_below_80 = mean(z < 80),
        p_e0_at_least_80 = mean(z >= 80),
        posterior_median = median(z),
        q025 = unname(quantile(z, 0.025, type = 8)),
        q975 = unname(quantile(z, 0.975, type = 8))
      )
    }
  }
}

state_results <- bind_rows(probability_rows) %>%
  mutate(
    health_index_group = factor(
      health_index_group,
      levels = c("Front-runner", "Achiever", "Aspirant")
    ),
    state = factor(state, levels = state_order),
    sex = factor(sex, levels = c("Female", "Male"))
  ) %>%
  arrange(health_index_group, state, sex, year)

# Long format for plotting and group-level summaries.
heatmap_data <- state_results %>%
  mutate(
    state = factor(as.character(state), levels = rev(state_order)),
    year_f = factor(as.character(year), levels = c("2024", "2050")),
    probability = p_e0_below_80,
    label = scales::percent(probability, accuracy = 0.1),
    label_colour = if_else(probability < 0.35, "white", "black")
  )

# -----------------------------------------------------------------------------
# 5. Heatmap: one threshold, faceted by sex; 2024 and 2050 side by side
# -----------------------------------------------------------------------------

n_aspirant <- sum(state_groups == "Aspirant")
n_achiever <- sum(state_groups == "Achiever")
separator_y <- c(n_aspirant + 0.5, n_aspirant + n_achiever + 0.5)

p_heatmap <- ggplot(heatmap_data, aes(x = year_f, y = state, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.35, width = 0.92, height = 0.92) +
  geom_text(aes(label = label, colour = label_colour), size = 2.35, show.legend = FALSE) +
  geom_hline(yintercept = separator_y, colour = "grey35", linewidth = 0.45) +
  facet_grid(. ~ sex) +
  scale_x_discrete(drop = FALSE, expand = expansion(add = 0.02)) +
  scale_y_discrete(drop = FALSE, expand = expansion(add = 0.02)) +
  scale_fill_viridis_c(
    option = "A",
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    name = "Posterior\nprobability"
  ) +
  scale_colour_identity() +
  labs(
    title = "Posterior probability of life expectancy below 80 years, 2024 and 2050",
    subtitle = "State-specific latent trajectories; states ordered by Health Index group",
    x = NULL,
    y = NULL,
    caption = paste(
      "Each tile shows P(e₀ < 80); P(e₀ ≥ 80) is its complement (100% minus the tile value).",
      "For states without a 2024 observation, the latent state is propagated from the latest observed year.",
      "2024 and 2050 probabilities use a common 0–100% colour scale."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(face = "bold", colour = "grey20"),
    axis.text.y = element_text(size = 7.2, colour = "grey15"),
    panel.spacing.y = grid::unit(0.5, "lines"),
    panel.spacing.x = grid::unit(0.8, "lines"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(size = 7.5, hjust = 0)
  )

ggsave(
  filename = file.path(figure_dir, "probability_below_80_heatmap_2024_2050.png"),
  plot = p_heatmap,
  width = 11,
  height = 12,
  units = "in",
  dpi = 400,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 6. Summary across Health Index groups
# -----------------------------------------------------------------------------
# These are equal-weighted summaries of state-specific probabilities, not
# population-weighted probabilities for the group population. The probability
# of e0 >= 80 is the complement of the probability of e0 < 80.

group_summary <- heatmap_data %>%
  group_by(year_f, sex, health_index_group) %>%
  summarise(
    n_states = n_distinct(state),
    mean_state_probability_below_80 = mean(probability),
    median_state_probability_below_80 = median(probability),
    min_state_probability_below_80 = min(probability),
    max_state_probability_below_80 = max(probability),
    n_states_probability_below_80_over_50 = sum(probability > 0.5),
    pct_states_probability_below_80_over_50 = 100 * mean(probability > 0.5),
    .groups = "drop"
  ) %>%
  mutate(
    year = as.integer(as.character(year_f)),
    mean_state_probability_at_least_80 = 1 - mean_state_probability_below_80,
    median_state_probability_at_least_80 = 1 - median_state_probability_below_80
  ) %>%
  select(-year_f) %>%
  relocate(year, sex, health_index_group)

write_csv(
  group_summary %>% mutate(health_index_group = as.character(health_index_group), sex = as.character(sex)),
  file.path(table_dir, "health_index_group_probability_below_80_2024_2050.csv")
)

# Save a pared-down state-level table with the two complementary outcomes.
state_probability_output <- state_results %>%
  transmute(
    health_index_group = as.character(health_index_group),
    state = as.character(state),
    sex = as.character(sex),
    year,
    p_e0_below_80,
    p_e0_at_least_80,
    posterior_median,
    q025,
    q975
  )

write_csv(
  state_probability_output,
  file.path(table_dir, "state_probability_below_80_2024_2050.csv")
)

cat("\nHEATMAP CREATED\n")
cat("Figure: ", file.path(figure_dir, "probability_below_80_heatmap_2024_2050.png"), "\n", sep = "")
cat("State-level probabilities: ", file.path(table_dir, "state_probability_below_80_2024_2050.csv"), "\n", sep = "")
cat("Health Index summaries: ", file.path(table_dir, "health_index_group_probability_below_80_2024_2050.csv"), "\n", sep = "")
cat("\n2024 estimates for states without a 2024 observation are model-based extrapolations from the state's latest observation.\n")
